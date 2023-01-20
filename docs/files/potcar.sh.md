```bash
#!/bin/bash    
# Create a GGA_PAW POTCAR file by concatenation of POTCAR files   
# Define local potpaw_GGA pseudopotential repository:
repo="<PATH TO potpaw_PBE FOLDER>"

# Check if older version of POTCAR is present
if [ -f POTCAR ] ; then                                                                   
	mv -f POTCAR old-POTCAR
	echo " ** Warning: old POTCAR file found and renamed to 'old-POTCAR'."
fi

# Main loop - concatenate the appropriate POTCARs (or archives)
for i in $*
do
	if test -f $repo/$i/POTCAR ; then
		cat $repo/$i/POTCAR >> POTCAR
	elif test -f $repo/$i/POTCAR.Z ; then   
		zcat $repo/$i/POTCAR >> POTCAR   
	elif test -f $repo/$i/POTCAR.gz ; then
		gunzip -c $repo/$i/POTCAR.gz >> POTCAR
	else
		echo " ** Warning: No suitable POTCAR for element '$i' found!! Skipped this element."
	fi       
done
```