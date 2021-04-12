#!/bin/bash
FSLPARALLEL=36;export FSLPARALLEL

######################## PROCESSING START ##########################################################################################
echo ""
echo " ****************************************************************"
echo " ---- group level analysis: randomise `date '+%Y-%m-%d %H:%M'`"
echo " ----------------------------------------------------------------"

iperm=5000
igroup=data
#iseed=mPFC
BOLD_RUN=func_resting2
processing_type=smooth5mm_denoised_mni_epirnorm
design_type=design_age

projectDir=/Volumes/andlab/project_AU
mainoutDir=${projectDir}/analysis; if [ ! -e $mainoutDir ]; then mkdir -p ${mainoutDir}; fi
designDir=${projectDir}/scripts/group_design/${design_type}
dataDir=${projectDir}/subs/${igroup}
covariateDocu=${projectDir}/documents/Andlab_ABIDE_basic_subject_Info_only_available_data.csv

maskDir=${projectDir}/mask
imask=${maskDir}/all_group_mask.nii.gz
standard_T1_mask=/Volumes/andlab/resources/brain_templates/standard/adults/MNI152_T1_2mm_brain_mask.nii.gz

for iseed in mPFC LC salience; do

	finaloutDir=${mainoutDir}/fcseed/${BOLD_RUN}/${processing_type}/${iseed}/groupLevel_wholebrain/randomise_${design_type}
	finalIndex=${finaloutDir}/done.group_processing

	if [ ! -d ${finaloutDir} ] || [ ! -e ${finalIndex} ]; then 
		mkdir -p ${finaloutDir}
		echo " ---- prep files for randomise"
		ls -f ${dataDir}/AUC*/analysis/fcseed/${BOLD_RUN}/${processing_type}/${iseed}/fcmap_zvalue.nii.gz > ${finaloutDir}/sublist.txt
		ls -f ${dataDir}/AUP*/analysis/fcseed/${BOLD_RUN}/${processing_type}/${iseed}/fcmap_zvalue.nii.gz >> ${finaloutDir}/sublist.txt
		sublist=${finaloutDir}/sublist.txt
		fslmerge -a ${finaloutDir}/4d_zmaps.nii.gz `cat ${finaloutDir}/sublist.txt`

		## generate group mask
		if [ ! -e ${imask} ]; then
			echo " ---- mask creation"
			echo "      : individual mni mask creation"
			cd ${dataDir}     
		    for isub in AU*; do
		    	if [ ! -e ${dataDir}/${isub}/raw/${BOLD_RUN}/${BOLD_RUN}_mni.mask.nii.gz ]; then
			    	inativemask=${dataDir}/${isub}/raw/${BOLD_RUN}/${BOLD_RUN}_mask.nii.gz
			    	iregDir=${dataDir}/${isub}/raw/${BOLD_RUN}/reg
			    	imat1=${iregDir}/ants.example_func2standard0GenericAffine.mat
					iwarp=${iregDir}/ants.example_func2standard1Warp.nii.gz
					istandard=${iregDir}/standard.nii.gz
					itransMats="-t ${imat1} -t ${iwarp}"

					antsApplyTransforms -d 3 -r ${istandard} -i ${inativemask} ${itransMats} \
					-o ${dataDir}/${isub}/raw/${BOLD_RUN}/${BOLD_RUN}_mni_mask.nii.gz
					
					fslmaths ${dataDir}/${isub}/raw/${BOLD_RUN}/${BOLD_RUN}_mni_mask.nii.gz \
					-bin ${dataDir}/${isub}/raw/${BOLD_RUN}/${BOLD_RUN}_mni_mask.nii.gz
				fi
			done

			echo "     : group mni mask creation"     
		    ## 2b. Merge the 2a brains using fslmerge -t
			fslmerge -t ${imask} `ls ${dataDir}/AU*/raw/${BOLD_RUN}/${BOLD_RUN}_mni_mask.nii.gz`
			
			## 2c. Make a group mask that contains only voxels present in all subjs using fslmaths -Tmin
			fslmaths ${imask} -Tmin ${imask};
			fslmaths ${imask} -bin ${imask}; 
			fslmaths ${imask} -mul ${standard_T1_mask} ${imask}; # exclude white matter
			cp ${imask} ${finaloutDir}
		else 
			cp ${imask} ${finaloutDir}
		fi

		## run randomise
		cd ${finaloutDir}
		echo " ---- run randomise TFCE (n =  ${iperm})"
		istart=`date '+%Y-%m-%d %H:%M'`
		randomise_parallel -i ${finaloutDir}/4d_zmaps.nii.gz -o ${finaloutDir}/TwoSampleT_TFCE \
		-d ${designDir}/design.mat -t ${designDir}/design.con -e ${designDir}/design.grp -m ${finaloutDir}/all_group_mask.nii.gz -n ${iperm} -T
		echo "TFCE start: ${istart} - done: `date '+%Y-%m-%d %H:%M'`" > ${finaloutDir}/done.group_processing
		echo " ****************************************************************"
		echo " ---- processing done: randomise TFCE `date '+%Y-%m-%d %H:%M'`"
		echo " ----------------------------------------------------------------"
		rm -rf *SEED*
		
		echo " ---- run randomise Cluster (n =  ${iperm})"
		istart=`date '+%Y-%m-%d %H:%M'`
		randomise_parallel -i ${finaloutDir}/4d_zmaps.nii.gz -o ${finaloutDir}/TwoSampleT_Cluster \
		-d ${designDir}/design.mat -t ${designDir}/design.con -e ${designDir}/design.grp -m ${finaloutDir}/all_group_mask.nii.gz -n ${iperm} -C 2.57
		echo "Cluster start: ${istart} - done: `date '+%Y-%m-%d %H:%M'`" >> ${finaloutDir}/done.group_processing
		echo " ****************************************************************"
		echo " ---- processing done: randomise Cluster `date '+%Y-%m-%d %H:%M'`"
		echo " ----------------------------------------------------------------"
		rm -rf *SEED*

		# echo " ---- run randomise Vox (n =  ${iperm})"
		# istart=`date '+%Y-%m-%d %H:%M'`
		# randomise_parallel -i ${finaloutDir}/4d_zmaps.nii.gz -o ${finaloutDir}/TwoSampleT_Vox \
		# -d ${designDir}/design.mat -t ${designDir}/design.con -m ${finaloutDir}/all_group_mask.nii.gz -n ${iperm} -x 
		# echo "Voxel start: ${istart} - done: `date '+%Y-%m-%d %H:%M'`" >> ${finaloutDir}/done.group_processing 
		# echo " ****************************************************************"
		# echo " ---- processing done: randomise Vox `date '+%Y-%m-%d %H:%M'`"
		# echo " ----------------------------------------------------------------"

		# rm -rf *SEED*
	fi
done
