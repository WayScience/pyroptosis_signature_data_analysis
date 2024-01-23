macro "Maximum Intensity Projection" {
	//INPUT/OUPUT folders
	//inDir=getDirectory("Choose the input folder");
	//outDir=getDirectory("Choose the ouput folder");
	inDir="/home/lippincm/Documents/ML/Interstellar_Analysis/8.cytopick_analysis/figures/PBMC/"
	outDir="/home/lippincm/Documents/ML/Interstellar_Analysis/8.cytopick_analysis/figures/PBMC/"
	myList=getFileList(inDir);  //an array
	start = getTime();
	waitForUser("I solemnly swear I am up to no good");
	// Make an array of tif files only
	flist = newArray(0);
	for (i=0; i<myList.length; i++) {
		if (endsWith(myList[i], ".tiff")) {
			flist = append(flist, myList[i]);
		}
	}

	for (j = 0 ; j < flist.length ; j++ ){
		progress = j / flist.length;
		progress = progress * 100;
		print(progress+"% complete");
		path=inDir+flist[j];
		open(path);
		a = getTitle();

		a = replace(a, ".tiff", "");

		run("Split Channels");
		selectImage("C1-"+a+".tiff");
		run("Magenta");
		selectImage("C2-"+a+".tiff");
		selectImage("C3-"+a+".tiff");
		run("Merge Channels...", "c2=C2-"+a+".tiff c5=C3-"+a+".tiff c6=C1-"+a+".tiff create");
		rename(a);
		saveAs("Tiff", outDir+a+".tif");
		saveAs("PNG", outDir+a+".png");
		close();

	}
	sec = (getTime()-start)/1000;
	min = sec/60;
	hour = min/60;
	print(sec+" seconds");
	print(min+" minutes");
	print(hour+" hours");
	waitForUser("Mischeif Managed");
}


function append(arr, value) {
    arr2 = newArray(arr.length+1);
    for (i=0; i<arr.length; i++)
        arr2[i] = arr[i];
        arr2[arr.length] = value;
    return arr2;
}
