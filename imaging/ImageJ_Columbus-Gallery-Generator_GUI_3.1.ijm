// -------------------------------------------------------------------
// Written by: Patrice Mascalchi
// Date: 2015
// Location: Cancer Research Institute, University of Cambridge, UK
// Contact: patrice.mascalchi@gmail.com
// -------------------------------------------------------------------

// version 2.0/2.1: automation to create galleries for each well
// 2.2: better adjustment of Min / Max for the display
// 2.5: Drawing circles around spots with choice of the channel
// to deactivate drawing: circled > nbch...
// 2.5b: !!!! no crop of the canvas !!! / same scale for intensities by default
// 3.0: generate merged images from all the fields in 1 well and draws circles + field number on it
// 3.1: create a new channel for annotations (to be able to disable the overlay)

macro "gallery-generator" {
// radius of circle in px (0 = no circle)
rad = 15;
radcol = "Cyan";

// 0 = no label
label = 1;

// nb of fields
nbfi = 28;

// channels colors
keepch = newArray(0, 1, 1, 1, 1);
chchoice = newArray("Grays", "Green", "Yellow", "Red");
chminmax = newArray(0, 0, 250, 5, 220, 0, 220, 0, 220); 

imgpath="C:\\Users\\mascal01\\Documents\\UsersData\\Lovorka\\echo-test-1\\Meas01\\";
//imgpath="C:\\Users\\Public\\Documents\\SCREEN-1-130314\\";
//imgpath="C:\\Users\\Public\\Documents\\SCREEN-2-020714\\";
if (imgpath=="") imgpath = getDirectory("Choose the MEAS folder containing all the flex files for 1 plate");
if (imgpath=="") exit;

respath = File.openDialog("Choose the CSV file that contains the individual results");
if (respath=="") exit;
// ****************************************************************************************************************************
run("Close All");

// GENERATE IMAGES ////////////////////////////////////////////////////////
// output
if (indexOf(imgpath, "[")>0 || indexOf(imgpath, "]")>0) exit("error with directory - replace \[ or \] from the folder's names");
outfld = imgpath + "field-images\\";
if (!File.exists(outfld)) File.makeDirectory(outfld);

indwell = indexOf(respath, ".result.");
well = substring(respath, indwell+8, indwell+11);
well = replace(well, "\\[", "");
letter = toUpperCase(substring(well, 0, 1));
number = substring(well,1);
wellname = letter+number;

alphab = newArray("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P");
for (a=0;a<alphab.length;a++){
	if (letter==alphab[a]) {
		letconv = a+1;
		a = alphab.length;
	}		
}
// FLEX format = rrrcccfff.flex
i1 = fillzeros(letconv,3);
i2 = fillzeros(number,3);

if (!File.exists(outfld + wellname + "_010.tif")) {
	setBatchMode(true);
	// loop *********************************************************************************
	run("Bio-Formats Macro Extensions");
	for (fi=1;fi<=nbfi;fi++) {
		run("Close All");
		i3 = fillzeros(fi,3);
		imgname = i1 + i2 + i3 + ".flex";
	
		//open(imgpath + imgname);
		run("Bio-Formats Importer", "open=["+ imgpath+imgname +"] autoscale color_mode=Default concatenate_series open_files open_all_series view=Hyperstack stack_order=Default");
		run("Make Composite", "display=Composite");
		run("8-bit");
		for (c=chchoice.length; c>0; c--) {
			setSlice(c);
			if (keepch[c]==0) {
				run("Delete Slice", "delete=channel");
			} else {
				setMetadata("Label", fi);
				setMinAndMax(chminmax[c*2-1], chminmax[c*2]);
				run("Apply LUT");
				//run(chchoice[c-1]);
			}
		}
		setSlice(nSlices);
		if (rad>0) run("Add Slice", "add=channel");

		//run("Stack to RGB");
		// Labelling
		if (label) {
			setFont("SanSerif", 10);
			makeText(fi, 5, 5);
			run("Draw", "stack");
		}
		saveAs("Tiff", ""+outfld+wellname+"_"+fillzeros(fi,2)+".tif");
		run("Close All");
	}
} // end of IF

closeWin("Log");
run("Close All");
setBatchMode(false);
 ///////////////////////////////////////////////////////////////////////////

resdir = File.getParent(respath)+"\\";
resname = substring(respath, indwell+19);

// output
if (indexOf(resdir, "[")>0 || indexOf(resdir, "]")>0) exit("error with directory - replace \[ or \] from the folder's names");
outdir2 = resdir + "galleries\\";
if (!File.exists(outdir2)) File.makeDirectory(outdir2);

data = File.openAsString(respath);
lines = split(data, "\n");
cols = split(lines[0], ",");		//titles
cols[0] = "";				//to make 1 empty entry = no filtering

moreless = newArray("more than", "less than");
cho = newArray(16);			//choices

Dialog.create("Choose the parameters to filter your data");
Dialog.addMessage("______________________ all the parameters below are optional ____________________");
Dialog.addChoice("Parameter 1:", cols);
Dialog.addChoice("should be ", moreless);
Dialog.addNumber("",0);
//Dialog.addMessage("______________________");
Dialog.addChoice("Parameter 2:", cols);
Dialog.addChoice("should be ", moreless);
Dialog.addNumber("",0);
//Dialog.addMessage("______________________");
Dialog.addChoice("Parameter 3:", cols);
Dialog.addChoice("should be ", moreless);
Dialog.addNumber("",0);
//Dialog.addMessage("______________________");
Dialog.addChoice("Parameter 4:", cols);
Dialog.addChoice("should be ", moreless);
Dialog.addNumber("",0);
//Dialog.addMessage("______________________");
Dialog.show;

setBatchMode(true);

for (i=0;i<4;i++) {
	ind = 4*i;
	cho[ind] = Dialog.getChoice;
	cho[ind+1] = Dialog.getChoice;
	cho[ind+2] = Dialog.getNumber;
}

// finding the columns numbers for coordinates and field image position
for (c=0;c<cols.length;c++) {
	if (cols[c]=="X") xC = c;
	if (cols[c]=="Y") yC = c;
	if (cols[c]=="Row") rowC = c;
	if (cols[c]=="Column") colC = c;
	if (cols[c]=="Field") fieldC = c;
}

data = File.openAsString(respath);
lines = split(data, "\n");
cols1 = split(lines[1], ",");

if (lines.length < 2) exit("No results in that table");

run("Bio-Formats Macro Extensions");

i1 = fillzeros(cols1[rowC],3);
i2 = fillzeros(cols1[colC],3);

spX = newArray(200);		// limited to 200 OBJECTS per image
spY = newArray(200);
objC = 0;
nSl = keepch[1]+keepch[2]+keepch[3]+keepch[4];
ifrad = 0;

for (i=1;i<lines.length;i++) {
	run("Close All");

	cls = split(lines[i], ",");
	if (i<lines.length-1) {
		clsnext = split(lines[i+1], ",");
	} else {
		clsnext = newArray(cols1.length);
	}
	
	// FLEX format = rrrcccfff.flex
	i3 = fillzeros(cls[fieldC],3);
	imgname = i1 + i2 + i3 + ".flex";

	if (parseInt(clsnext[fieldC])==parseInt(cls[fieldC]) && i!=lines.length-1) {	// do not load the image yet
		objC = objC+1;
		spX[objC] = cls[xC];
		spY[objC] = cls[yC];
	} else {
		open(outfld+wellname+"_"+i3+".tif");			
		// CIRCLES DRAWING --------
		if (rad>0) {
			setSlice(nSl+1);
			setForegroundColor(255, 255, 255);
		    	spX[objC+1] = cls[xC];
		    	spY[objC+1] = cls[yC];
		    	//print(c+"--"+spX[objC]+"/"+spY[objC]);		    	
			for (sp=1;sp<=objC+1;sp++) {
				makeOval(spX[sp]-rad, spY[sp]-rad, rad*2, rad*2);
				run("Draw", "slice");
				run("Select None");
			}
			ifrad = 1;
		}
		objC=0;		//reset
		saveAs("Tiff", ""+outfld+wellname+"_"+i3+".tif");
		run("Close All");
	}
}	// end for

// generate stack
run("Image Sequence...", "open=["+outfld+wellname+"_001.tif] file="+wellname+" sort");
run("Stack to Hyperstack...", "order=xyczt(default) channels="+nSl+ifrad+" slices=1 frames="+nbfi+" display=Composite");
for (c=1;c<=nSl;c++) {
	setSlice(c);
	run(chchoice[c-1]);
}
setSlice(nSl+ifrad);
run(radcol);

saveAs("Tiff", ""+outdir2+wellname+"_"+resname+".tif");
wait(100);
for (fi=1;fi<=nbfi;fi++) File.delete(outfld+wellname+"_"+fillzeros(fi,2)+".tif");
closeWin("Log");

run("Close All");

setBatchMode(false);
open(""+outdir2+wellname+"_"+resname+".tif");
setLocation(100, 100);
run("Channels Tool...");
selectWindow("Channels");
setLocation(800, 100);
run("Brightness/Contrast...");
selectWindow("B&C");
setLocation(1000, 100);
Stack.setActiveChannels("10011");

} // end of macro --------------------------------------------------------------------------------------------

}//-----------------------------------------------
function fillzeros(nb, zeronb) {
	for (z=0;z<=zeronb;z++) nb = "0" + nb;
	n = substring(nb, lengthOf(nb)-3);
	return n;
}//-----------------------------------------------
function closeWin(name) {
	if (isOpen(name)) {
	     selectWindow(name);
	     run("Close");
	}
}