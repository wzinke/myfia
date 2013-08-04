MyFIA
=====

Miscellaneous bash scripts to facilitate the processing and handling of fMRI data. The scripts primarily utilize tools provided by the FSL package (www.fmrib.ox.ac.uk/fsl).

I generated these scripts while working with primate fMRI data and collected them as my *MaFIA* (Macaque Functional Image Analysis) toolbox. But most of the scripts and tools do not care about the species the data is coming from, thus several scripts might be just as useful for human fMRI data. Therefore, I selected several scripts and included them in my funcional image analysis package (*MyFIA*). I hope they are helpful for other users. However, I can not guarantee that they are fully functional and free of bugs.

To use the scripts, download single files or clone the complete directory:

    git clone https://github.com/wzinke/myfia.git

and then add the directory to the known PATHs, i.e. in the ~/.bash_profile file:

    PATH=<path_to_dir>/myfia:${PATH}
    export  PATH

Most of the BASH scripts require at least the FSL package (www.fmrib.ox.ac.uk/fsl), some need AFNI (http://afni.nimh.nih.gov/afni/).

The scripts are released under the MIT licence, see MIT-LICENSE.txt for details. I do not guarantee that they are free of bugs, on the contrary, I am happy for any bug reports or comments of possible improvements.
