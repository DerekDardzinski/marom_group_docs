# potcar.sh
The `potcar.py` file can be used to generate kpoints for all possible types of calculations that we do in the Marom group. To use it, copy the code block below and put it in `~/bin/potcar.py` then run:

```bash
chmod +x ~/bin/potcar.py
```

Additionally, you will need to add the path to your `potpaw_PBE` file which holds the VASP pseudopotential information.

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