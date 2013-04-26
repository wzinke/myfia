MyFIA
=====

Miscellaneous bash scripts to facilitate the processing and handling of fMRI data, primarily based on tools provided by the FSL package.
I generated these scripts while working with primate fMRI data and collected them in my MaFIA (Macaque Functional Image Analysis). But most of the times, the scripts and tools do not care about the species the data is coming from, thus several scripts might be just as useful for human fMRI data. Therefore, I selected several of my scripts and included them in my funcional image analysis package (MyFIA). I hope they are helpful for other users. However, I can not guarantee that they are fully functional and free of bugs.

To use the scripts, add the directory to the known PATHs, i.e. in the ~/.bash_profile file:

    PATH=<path_to_dir>/myfia:${PATH}
    export  PATH

Most of the BASH scripts require the FSL package (www.fmrib.ox.ac.uk/fsl).

The scripts are released under the MIT licence, see MIT-LICENSE.txt for details.
