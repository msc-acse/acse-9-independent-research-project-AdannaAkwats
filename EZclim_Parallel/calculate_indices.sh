#!/bin/bash
func=$1
orig=$2
base=$3
timely=$4
save=$5 
res=$6

echo "CDO Script running ---------------------------------------------------------------------"
echo "Arguments passed:"
echo " - index: "$func 
echo " - future run: " $orig
echo " - control run: " $base
echo " - time step: "$timely
echo " - output file: "$save
echo " - output directory: "$res

# Calculate ENSO index 
if [ $func = "enso" ]
then
	if [ $timely = "daily" ]; then
  		cdo ydaysub -fldmean -sellonlatbox,-170,-120,-5,5 $orig -ydaymean -fldmean -sellonlatbox,-170,-120,-5,5 $base SSTAnom.nc
	else
  		cdo ymonsub -fldmean -sellonlatbox,-170,-120,-5,5 $orig -ymonmean -fldmean -sellonlatbox,-170,-120,-5,5 $base SSTAnom.nc
	fi
	cdo div SSTAnom.nc -timstd SSTAnom.nc SSTAnomStd.nc
	cdo runmean,5 SSTAnomStd.nc Nino34.nc
	cdo runmean,3 SSTAnomStd.nc $save
# Delete all files except final output
rm SSTAnom.nc
rm SSTAnomStd.nc
rm Nino34.nc
fi	

# Calculate TNI index and NINO 1+2 index
if [ $func = "nino12" ] || [ $func = "tni" ]
then
nino12="NINO12.nc"
resnino12="${res}${nino12}"
	if [ $timely = "daily" ]; then
		cdo ydaysub -fldmean -sellonlatbox,-90,-80,-10,0 $orig -ydaymean -fldmean -sellonlatbox,-90,-80,-10,0 $base $resnino12
	else
		cdo ymonsub -fldmean -sellonlatbox,-90,-80,-10,0 $orig -ymonmean -fldmean -sellonlatbox,-90,-80,-10,0 $base $resnino12
	fi
fi

# Calculate TNI index and NINO 4 index 
if [ $func = "nino4" ] || [ $func = "tni" ]
then
nino4="NINO4.nc"
resnino4="${res}${nino4}"
$
	if [ $timely = "daily" ]; then 
		cdo ydaysub -fldmean -sellonlatbox,160,-150,-5,5 $orig -ydaymean -fldmean -sellonlatbox,160,-150,-5,5 $base $resnino4
	else
		cdo ymonsub -fldmean -sellonlatbox,160,-150,-5,5 $orig -ymonmean -fldmean -sellonlatbox,160,-150,-5,5 $base $resnino4
	fi
fi

# Use NINO 1+2 and NINO 4 to calculate TNI index
if [ $func = "tni" ]; then
	cdo div $resnino12 -timstd $resnino12 Nino12AnomStd.nc
	cdo div $resnino4 -timstd $resnino4 Nino4AnomStd.nc
	cdo sub Nino12AnomStd.nc Nino4AnomStd.nc TNIRaw.nc
	cdo runmean,5 TNIRaw.nc TNISmth.nc
	cdo div TNISmth.nc -timstd TNISmth.nc $save

# Remove all files except final output
rm $resnino4
rm Nino4AnomStd.nc
rm $resnino12
rm Nino12AnomStd.nc
rm TNIRaw.nc
rm TNISmth.nc
fi

# Calculate IOD index
if [ $func = "iod" ]; then
	if [ $timely = "daily" ]; then 
		cdo ydaysub $orig -ydaymean $base SSTAnom.nc
	else
		cdo ymonsub $orig -ymonmean $base SSTAnom.nc
	fi
	cdo sub -fldmean -sellonlatbox,50,70,-10,10 SSTAnom.nc -fldmean -sellonlatbox,90,110,-10,0 SSTAnom.nc $save
rm SSTAnom.nc
fi

# Calculate AMO index
if [ $func = "amo" ]; then
        if [ $timely = "daily" ]; then
        	cdo ydaysub -sellonlatbox,-180,180,-60,60 $orig -ydaymean -sellonlatbox,-180,180,-60,60 $base GlbAnom.nc
		cdo ydaysub -sellonlatbox,-120,0,0,60 $orig -ydaymean -sellonlatbox,-120,0,0,60 $base AtlAnom.nc
        else
		cdo ymonsub -sellonlatbox,-180,180,-60,60 $orig -ymonmean -sellonlatbox,-180,180,-60,60 $base GlbAnom.nc
		cdo ymonsub -sellonlatbox,-120,0,0,60 $orig -ymonmean -sellonlatbox,-120,0,0,60 $base AtlAnom.nc
        fi
        cdo runmean,10 -yearmean -fldmean GlbAnom.nc GlbbAnnSmth.nc
	cdo runmean,10 -yearmean -fldmean AtlAnom.nc AtlAnnSmth.nc
	cdo sub AtlAnnSmth.nc GlbbAnnSmth.nc $save
rm GlbAnom.nc
rm AtlAnom.nc
rm GlbbAnnSmth.nc
rm AtlAnnSmth.nc
fi

# Calculate PDO index
if [ $func = "pdo" ]; then
	if [ $timely = "daily" ]; then
		cdo ydaysub -sellonlatbox,110,-100,20,70 $orig -ydaymean -sellonlatbox,110,-100,20,70 $base SSTAnom.nc
	else
		cdo ymonsub -sellonlatbox,110,-100,20,70 $orig -ymonmean -sellonlatbox,110,-100,20,70 $base SSTAnom.nc
	fi
	cdo detrend SSTAnom.nc SSTAnomDtrnd.nc
	cdo fldmean -sellonlatbox,-180,180,-60,70 SSTAnomDtrnd.nc GlbMeans.nc
	cdo sub SSTAnomDtrnd.nc -enlarge,SSTAnomDtrnd.nc GlbMeans.nc ResAnom.nc
	cdo eof,1 ResAnom.nc eigenvalues.nc eigenvectors.nc
	cdo eofcoeff eigenvectors.nc ResAnom.nc $save
rm SSTAnom.nc
rm SSTAnomDtrnd.nc
rm GlbMeans.nc
rm ResAnom.nc
rm eigenvalues.nc
rm eigenvectors.nc
fi

# Calculate AAO index
if [ $func = "aao" ]; then
	cdo -sellevel,700 $orig level700.nc
	cdo -sellevel,700 $base level700Base.nc
	if [ $timely = "daily" ]; then
		cdo ydaysub -sellonlatbox,-180,180,-90,-20 level700.nc -ydaymean -sellonlatbox,-180,180,-90,-20 level700Base.nc 700Anom.nc
	else
		cdo ymonsub -sellonlatbox,-180,180,-90,-20 level700.nc -ymonmean -sellonlatbox,-180,180,-90,-20 level700Base.nc 700Anom.nc
	fi
	cdo eof,1 700Anom.nc eigenvalues.nc eigenvectors.nc
	cdo eofcoeff eigenvectors.nc 700Anom.nc $save
rm level700.nc
rm level700Base.nc
rm 700Anom.nc
rm eigenvalues.nc
rm eigenvectors.nc
fi

# Calculate AO index
if [ $func = "ao" ]; then
	cdo -sellevel,1000 $orig level1000.nc
        cdo -sellevel,1000 $base level1000Base.nc
        if [ $timely = "daily" ]; then
                cdo ydaysub -sellonlatbox,-180,180,-90,-20 level1000.nc -ydaymean -sellonlatbox,-180,180,-90,-20 level1000Base.nc 1000Anom.nc
        else
                cdo ymonsub -sellonlatbox,-180,180,-90,-20 level1000.nc -ymonmean -sellonlatbox,-180,180,-90,-20 level1000Base.nc 1000Anom.nc
	fi
	cdo eof,1 1000Anom.nc eigenvalues.nc eigenvectors.nc
	cdo eofcoeff eigenvectors.nc 1000Anom.nc $save
rm level1000.nc
rm level1000Base.nc
rm 1000Anom.nc
rm eigenvectors.nc
rm eigenvalues.nc
fi

# Calculate NAO index
if [ $func = "nao" ]; then
	if [ $timely = "daily" ]; then
		cdo ydaysub -sellonlatbox,-90,40,20,70 $orig -ydaymean -sellonlatbox,-90,40,20,70 $base SLPAnom.nc
	else
		cdo ymonsub -sellonlatbox,-90,40,20,70 $orig -ymonmean -sellonlatbox,-90,40,20,70 $base SLPAnom.nc
	fi
	cdo eof,1 SLPAnom.nc eigenvalues.nc eigenvectors.nc
	cdo eofcoeff eigenvectors.nc SLPAnom.nc $save
rm SLPAnom.nc
rm eigenvalues.nc
rm eigenvectors.nc
fi
		
echo "done"
