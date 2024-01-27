macro "Maximum Intensity Projection" {
	//INPUT/OUPUT folders
	inDir=getDirectory("Choose the input folder");
	outDir=getDirectory("Choose the ouput folder");
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


		print(a);
		if (indexOf(a, "blue") >= 0) {
			selectWindow(a);
			run("Blue");
			print("match");

		} else if (indexOf(a, "green") >= 0) {
			selectWindow(a);
			run("Green");
			print("match");
		} else if (indexOf(a, "red") >= 0) {
			selectWindow(a);
			run("Red");
			print("match");
		} else if (indexOf(a, "yellow") >= 0) {
			selectWindow(a);
			run("Yellow");
			print("match");
		} else if (indexOf(a, "magenta") >= 0) {
			selectWindow(a);
			run("Magenta");
			print("match");
		} else {
			print("no matches found");
		}

		selectWindow(a);
		a = replace(a, ".tiff", "");
		saveAs("Tiff", outDir+a+".tif");
		saveAs("PNG", outDir+a+".png");
		selectWindow(a+".png");
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
