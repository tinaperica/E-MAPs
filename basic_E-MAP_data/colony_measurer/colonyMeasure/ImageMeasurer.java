/*
 * Created on Mar 6, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package colonyMeasure;

/**
 * @author Sean
 *
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
 
 // assen added tons of retarded info dumped on the console
 // todo: implement debug level static option to switch them on and off

import java.util.Arrays;
import java.util.Collections;

import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.measure.ResultsTable;
import ij.measure.ResultsTablePlus;
import ij.plugin.filter.ImageTools;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.RGBStackSplitterSean;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.plugin.filter.PixelMean;
import ij.io.FileSaver;

import java.awt.Rectangle;
import java.io.*;

public class ImageMeasurer
	{
	static int SHOW_NOTHING=0;
	static int SHOW_GRID_OVERLAY=1;

	// added by Assen, this would be the suffix of the grid overlay image
	private static final String MASKED_FILE_SUFFIX = "_grid_ovr.jpg";
	
	double PLATE_EDGE_HEIGHT_FRAC;
	double PLATE_EDGE_WIDTH_FRAC;
	int MINIMUM_PLATE_HEIGHT;
	int MINIMUM_PLATE_WIDTH;
	
	public int THRESHOLD = 0;
	
	//Image cropping details
	int cropBot,cropTop,cropL,cropR;

	ResultsTablePlus rtTotal,rtF,rt2;
	ColonyGrid grid, orig_grid;
	
	int ERROR_FLAG = ColMeasureProgram.ERR_NONE;			// default, everything is ok
	
	private MeasuredResults imageResults;
	
	private ImagePlus im, im_red1, im_redTemp;
	private ImageProcessor ip_red_ch;
	private String filename;
	private ImageMeasurerOptions imo;
	
	public MeasuredResults testResults;
	
	// constructor
	public ImageMeasurer(ImageProcessor ip, String fn, ImageMeasurerOptions im_o)
		{
		filename=fn;
		imo  = im_o;
		THRESHOLD = imo.threshold;
		
		ip_red_ch = ip;
		im = new ImagePlus("temp",ip);
		
		// extract the red chanel
		if (im.getBitDepth()==24)
			{ ip_red_ch = getRedChannel(ip); }
		else
			{ ip_red_ch = ip.convertToByte(false); }
		
		// initialize those two images
		im_red1 = new ImagePlus("Red channel",ip_red_ch);
		im_redTemp = new ImagePlus("Red channel",ip_red_ch);  	//	used to figure out how to rotate the image
		}
	
		
	public void setThreshold()
		{
		int initial_thres = THRESHOLD;
		if (initial_thres != 0) { System.out.println("ImageMeasurer->setThreshold: Threshold already set to " + initial_thres + " Not setting."); }	
		// mode 0 for de novo, 1 for use precalculated		
		setCroppingParameters(imo.gridType);
		im_redTemp = doImageCropping(imo, im_redTemp);
	
		if (initial_thres == 0)
			{
			System.out.println("ImageMeasurer->setThreshold: Threshold not set. Setting.");	 
			THRESHOLD = getColonyThreshold(im_redTemp); System.out.println("im_redTemp T: "+THRESHOLD);
			}
			
		rtTotal = getParticles(im_redTemp,THRESHOLD,50,5000);
		createColonyGrid(imo, im_redTemp);  //sets the value of grid
		im_redTemp.flush(); //attempting to minimize memory wastage
		
		// *** Rotate the image to try to account for the grid slant, and then set up the grid again
		im_red1 = new ImagePlus("Red channel", performEstimatedRotation(im_red1,grid.slope));
		im_redTemp = im_red1; 		 //redTemp will store the uncropped image in case we need it later
		im_red1 = doImageCropping(imo, im_red1);
		
		if (initial_thres == 0)
			{
			THRESHOLD = getColonyThreshold(im_red1);
			System.out.println("ImageMeasurer->setThreshold: im_red1 T: "+THRESHOLD);
			System.out.println("ImageMeasurer->setThreshold FINAL: " + THRESHOLD);
			}
					
		rtTotal=getParticles(im_red1,THRESHOLD,50,5000);
		createColonyGrid(imo, im_red1);  //sets the value of grid	
		}							


	public void fullAnalysis(int mode)
		{
		int j;
		
		// spit out some info
		System.out.println("ImageMeasurer->fullAnalysis mode: "+mode);
		System.out.println("ImageMeasurer->fullAnalysis output: " + filename);
		System.out.println("ImageMeasurer->fullAnalysis Threshold: " + THRESHOLD);
		System.out.println("Java memory in use = " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));

		setThreshold();
		
		// if mode is 0 just calc the threshold and return
		if (mode ==  ColMeasureProgram.AN_THRES_ONLY)
			{ System.out.println("No analysis, done on " + filename + ". Just setting threshold.");	 return; }
		
		// still zero, complain. unlikely to happen but still
		if (THRESHOLD == 0)
			{ System.out.println("ImageMeasurer->fullAnalysis: THRESHOLD is still 0. Dosen't sound right. Nothing done. Aborting."); return; }				
		
		// tell the world what the THRESHOLD now is
		System.out.println("ImageMeasurer->fullAnalysis: Using THRESHOLD " + THRESHOLD); 
		
		// *** Continue with analysis
		grid.applyGridMask(im_red1);
		rtF = getParticles(grid.maskedImage,THRESHOLD,10,5000);
		imageResults=new MeasuredResults(grid,rtF);
		orig_grid=grid;
		for (j=0;j<2;j++)
			{
			ERROR_FLAG=0;
			if ((imageResults.typicalNumSpots>1.1) | (imageResults.numzero>(0.33*grid.numCol*grid.numRow*grid.numDup)))
				{ //Take a second shot at it
				System.out.println("Taking a another try at the image");
				ERROR_FLAG=1; //if it gets here twice, it will report an error, even though it might still process the image correctly
				switch (imo.gridType)
					{
					case ColonyGrid.GRID_96: grid=new ColonyGrid96(im_red1,rtTotal,orig_grid,j); break;
					case ColonyGrid.GRID_384: grid=new ColonyGrid384(im_red1,rtTotal,orig_grid,j); break;
					case ColonyGrid.GRID_768a: grid=new ColonyGrid768a(im_red1,rtTotal,orig_grid,j); break;
					case ColonyGrid.GRID_1536: grid=new ColonyGrid1536(im_red1,rtTotal,orig_grid,j); break;
					}
				grid.applyGridMask(im_red1);
				rtF=getParticles(grid.maskedImage,THRESHOLD,10,5000);
				testResults=new MeasuredResults(grid,rtF);
				if (testResults.isBetterThan(imageResults)) imageResults=testResults;
				}
			}
	
		if (imageResults.gridExceedsImageBounds)
			{
			System.out.println("Original grid exceeded original cropped image - enlarging image and trying again");
			int numExtraPixels=(int)Math.round(im_red1.getHeight()*0.025);
			im_red1 = generateLargerCroppedImage(im_redTemp,numExtraPixels);
			grid = imageResults.grid;  //Retrieve the best grid we've had so far
			grid.adjustGridPositionForExpandedImage(numExtraPixels);
			grid.applyGridMask(im_red1);
			rtF=getParticles(grid.maskedImage,THRESHOLD,10,5000);
			imageResults = new MeasuredResults(grid,rtF);
			}
			
		imageResults.exportResults(filename);
		
		// error flagging
		if (imageResults.inconsistentReplicates()) { ERROR_FLAG = ColMeasureProgram.ERR_REPLICATES; }
			
		if ((rtF.numRows<(grid.numCol*grid.numRow*grid.numDup*0.75)) | (imageResults.numzero>(0.33*grid.numCol*grid.numRow*grid.numDup)))
				{ System.out.println("ImageMeasurer->errorFlagging: too few colonies"); ERROR_FLAG = ColMeasureProgram.ERR_TOO_FEW_COL; }
			
		if (((rtF.numRows>(grid.numCol*grid.numRow*grid.numDup*2)) | (imageResults.typicalNumSpots>1.1)))
			{ System.out.println("ImageMeasurer->errorFlagging: too many colonies"); ERROR_FLAG = ColMeasureProgram.ERR_TOO_MANY_COL; }
			
		if ((im_red1.getHeight()-grid.gridHeight>grid.heightTol))
				{ System.out.println("ImageMeasurer->errorFlagging: height tolerance error"); ERROR_FLAG = ColMeasureProgram.ERR_CROP_ERROR; }
			
		if (imageResults.gridExceedsImageBounds)
			{ System.out.println("ImageMeasurer->errorFlagging: cropping error - bounds of grid exceeded bounds of image"); ERROR_FLAG = ColMeasureProgram.ERR_CROP_ERROR; }
			
		if ((rtF.numRows < (grid.numCol*grid.numRow*grid.numDup*0.75)) |
			(rtF.numRows > (grid.numCol*grid.numRow*grid.numDup*2)) |
			(im_red1.getHeight()-grid.gridHeight>grid.heightTol) |
			(imageResults.typicalNumSpots>1.1) |
			imageResults.gridExceedsImageBounds |
			imageResults.inconsistentReplicates())
				{
				if (ERROR_FLAG == ColMeasureProgram.ERR_OTHER) 
					{ System.out.println("ImageMeasurer->errorFlagging: Misc Error processing file"); }
				}
		
		// some more setup	
		if (imo.showOption == SHOW_GRID_OVERLAY)
			{ imageResults.grid.maskedImage.show();	}
		
		// added by Assen
		if (imo.writeMaskedOption == 1)
			{	
			FileSaver fs = new FileSaver(imageResults.grid.maskedImage);
			String fname = filename+MASKED_FILE_SUFFIX;
			fs.saveAsJpeg(fname);
			System.out.println("ImageMeasurer->fullAnalysis: Writing grid overlay picture" + fname);
			}
		}



	private void createColonyGrid(ImageMeasurerOptions imo, ImagePlus im)
		{
		switch (imo.gridType)
			{
			case ColonyGrid.GRID_96: grid = new ColonyGrid96(im, rtTotal); break;
			case ColonyGrid.GRID_384: grid = new ColonyGrid384(im, rtTotal); break;
			case ColonyGrid.GRID_768a: grid = new ColonyGrid768a(im, rtTotal); break;
			case ColonyGrid.GRID_1536: grid = new ColonyGrid1536(im, rtTotal); break;
			}
		}


	private ImagePlus doImageCropping(ImageMeasurerOptions imo, ImagePlus im)
		{
		if ((imo.cropType==ColonyGrid.AUTO_CROP_OPTION)|(imo.cropType==ColonyGrid.AUTO_CROP_DECIDE_OPTION))
			{ im = new ImagePlus("Red channel",cropOutsidePlate(im,imo)); }
		else
			{ im = new ImagePlus("Red channel",justCropPlateEdges(im)); }
		return im;
		}
			
			
	protected static boolean spotIsNearCenter(ResultsTablePlus rt, int ind, OvalRoi roi)
		{
		int x=(int)Math.round(rt.getValue("X",ind));
		int y=(int)Math.round(rt.getValue("Y",ind));
		Rectangle r=roi.getBounds();
		double diameter=r.width;
		double radius=r.width*7/20; //the center of the colony should be inside a smaller oval
		double centX=r.x+diameter/2;
		double centY=r.y+diameter/2;
		boolean res = (distance(x,y,centX,centY)<radius);
		return res;
		}
	
	
	protected static double distance(double x1, double y1, double x2, double y2)
		{
		double distance=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
		distance=Math.sqrt(distance);
		return distance;
		}
			
	private void writeln(BufferedWriter bw, String str)
		{
		try { bw.write(str, 0, str.length()); bw.newLine(); }	
		catch (IOException io)
			{ System.out.println("Error writing a line to the file"); }
		}
		
		
	private ImageProcessor getRedChannel(ImageProcessor ip)
		{
		RGBStackSplitterSean splitter = new RGBStackSplitterSean();
		splitter.split(new ImagePlus("tempImage",ip));
		ImagePlus red=new ImagePlus("tempImage",splitter.red);
		return red.getProcessor();
		}

	
	protected static ResultsTablePlus getParticles(ImagePlus im, int thresh, int minSize, int maxSize)
		{
		ImagePlus im2=ImageTools.threshold(im,thresh);
		im2.setRoi(0,0,im2.getWidth(),im2.getHeight());
		ResultsTablePlus rt = new ResultsTablePlus();
		int options=0;
		int measurements=ParticleAnalyzer.AREA + ParticleAnalyzer.CENTROID+ParticleAnalyzer.ELLIPSE+ParticleAnalyzer.CIRCULARITY; //+ParticleAnalyzer.SHOW_OUTLINES;
		ParticleAnalyzer pa=new ParticleAnalyzer(options,measurements,rt,minSize,maxSize);
		pa.analyze(im2);
		rt.numRows=rt.getCounter();
		return rt;
		}
		
		
	protected static int[] findSpotsInOval(ResultsTablePlus rt, OvalRoi roi)
		{
		int list1[]=new int[rt.numRows];
		int x,y;
		int i,j=0;
		Rectangle r=roi.getBoundingRect();
		for (i=0;i<rt.numRows;i++)
			{
			x=(int)Math.round(rt.getValue("X",i));
			y=(int)Math.round(rt.getValue("Y",i));
			
			if (roi.contains(x,y))
				{ list1[j]=i; j=j+1; }
			}
			
		int list[]=new int[j];
		for (i=0;i<j;i++)
			{ list[i]=list1[i]; }

		return list;
		}
		
	
	public static int getColonyThreshold(ImagePlus im)
		{
		int h=im.getHeight();
		int w=im.getWidth();
		im.setRoi(w/4,h/4,w/2,h/2);
		ImagePlus im2=ImageTools.copy(im,im.getRoi());
		ImageProcessor ip=im2.getProcessor();

		return ((ByteProcessor)ip).getAutoThreshold();
		}

	
	protected static double getBaselinePixelValue(ImagePlus im)
		{
		int h=im.getHeight();
		int w=im.getWidth();
		double val1[] = new double[4];
		
		val1[0]=PixelMean.get(im,0,1,w,2);
		val1[1]=PixelMean.get(im,0,h-3,w,2);
		val1[2]=PixelMean.get(im,1,0,2,h);
		val1[3]=PixelMean.get(im,w-3,0,2,h);
		Arrays.sort(val1);
		
		return (val1[1]+val1[2])/2;
		}
		
		
	protected static int partialSum(int[] list, int start, int stop)
		{
		int i;
		int tot=0;
		if (start<0)
			{ start=0; }
			
		for (i=start;i<=stop;i++)
			{ tot=tot+list[i]; }

		return tot;
		}
		
		
	private ImageProcessor cropOutsidePlate(ImagePlus im, ImageMeasurerOptions imo)
		{
		double med = 0;
		int h=im.getHeight();
		int w=im.getWidth();
		int diff;

		double baseline=getBaselinePixelValue(im);
		cropBot = findVertCropPos(im, h-1, -1, baseline);
		cropTop = findVertCropPos(im, 0, 1, baseline);
		int tempFullHeight = cropBot-cropTop;
		cropBot =cropBot-(int)Math.round(PLATE_EDGE_HEIGHT_FRAC*tempFullHeight);
		cropTop = cropTop+(int)Math.round(PLATE_EDGE_HEIGHT_FRAC*tempFullHeight);
		diff = cropBot-cropTop;
		
		if ((diff < MINIMUM_PLATE_HEIGHT) & (imo.cropType == ColonyGrid.AUTO_CROP_DECIDE_OPTION))
			{
			System.out.println("ImageMeasurer->cropOutsidePlate: " + diff + " minimum plate height");
			cropTop=(int)Math.round(PLATE_EDGE_HEIGHT_FRAC*h);
			cropBot=h-(int)Math.round(PLATE_EDGE_HEIGHT_FRAC*h);
			}
			
		cropL = findHorCropPos(im,0,1, baseline);
		cropR = findHorCropPos(im,w-1,-1, baseline);
		int tempFullWidth=cropR-cropL;
		cropL = cropL + (int)Math.round(PLATE_EDGE_WIDTH_FRAC*tempFullWidth);
		cropR = cropR - (int)Math.round(PLATE_EDGE_WIDTH_FRAC*tempFullWidth);

		if ((cropR-cropL<MINIMUM_PLATE_WIDTH)&(imo.cropType==ColonyGrid.AUTO_CROP_DECIDE_OPTION))
			{
			cropL=(int)Math.round(w*PLATE_EDGE_WIDTH_FRAC);
			cropR=w-(int)Math.round(w*PLATE_EDGE_WIDTH_FRAC);
			}
			
		im.setRoi(cropL,cropTop,cropR-cropL,cropBot-cropTop);
		return im.getProcessor().crop();
		}
		
		
	private ImagePlus generateLargerCroppedImage(ImagePlus im, int numExtraPixels)
		{
		cropL=cropL-numExtraPixels;
		cropR=cropR+numExtraPixels;
		cropTop=cropTop-numExtraPixels;
		cropBot=cropBot+numExtraPixels;
		im.setRoi(cropL,cropTop,cropR-cropL,cropBot-cropTop);
		return new ImagePlus("Red channel",im.getProcessor().crop());
		}
		
		
	private ImageProcessor justCropPlateEdges(ImagePlus im)
		{
		double med=0;
		int h=im.getHeight();
		int w=im.getWidth();
		int cropBot,cropTop,cropL,cropR;

		cropBot = h;
		cropTop = 0;
		int tempFullHeight = cropBot-cropTop;
		cropBot=cropBot-(int)Math.round(PLATE_EDGE_HEIGHT_FRAC*tempFullHeight);
		cropTop=cropTop+(int)Math.round(PLATE_EDGE_HEIGHT_FRAC*tempFullHeight);
		cropL = 0;
		cropR = w;
		int tempFullWidth=cropR-cropL;
		cropL = cropL + (int)Math.round(PLATE_EDGE_WIDTH_FRAC*tempFullWidth);
		cropR = cropR - (int)Math.round(PLATE_EDGE_WIDTH_FRAC*tempFullWidth);
		im.setRoi(cropL,cropTop,cropR-cropL,cropBot-cropTop);
		return im.getProcessor().crop();
		}
		
		
	private ImageProcessor performEstimatedRotation(ImagePlus im, double slope)
		{
		ImageProcessor res;
		int h=im.getHeight();
		int w=im.getWidth();
		int cropTop=(int)Math.round(h*0.025);
		int cropH=(int)Math.round(h*0.95);
		int cropL=(int)Math.round(w*0.025);
		int cropW=(int)Math.round(w*0.95);
		double angle=Math.atan(slope)/Math.PI*-180;
		
		res=im.getProcessor();
		res.rotate(angle);
		res.setRoi(cropL, cropTop, cropW, cropH);
		System.out.println("ImageMeasurer->performEstimatedRotation: Rotating image");
		return res.crop();
		}
		
		
	protected static int findVertCropPos(ImagePlus im, int start,int dir, double baseline)
		{
		double med;
		int w=im.getWidth();
		int h=im.getHeight();
		int i;
		int[] medValues;
		int top=h/2 - 10;
		int[] tops;
		double typical=PixelMean.get(im,0,top,w,top+20);
		double thresh=(typical+baseline)/2;
		//System.out.println("baseline: "+baseline+" "+typical+" "+thresh);
		medValues=new int[Math.round(h/3)+1];
		tops=new int[Math.round(h/3)+1];
		for (i=0;i<h/4;i++)
			{
			top=start+i*1*dir+(dir-1);
			tops[i]=top;
			med=PixelMean.get(im,0,top,w,2);
			
			if (med>thresh)
				{ medValues[i]=1; }
			else
				{ medValues[i]=0; }
			}

		i=i-2;
		int flag=1;
		while (flag>0)
			{
			if (partialSum(medValues,i-9,i) == 0)
				{ flag=0; }
			i=i-1;
			}
			
		i=tops[i];
		//System.out.println("Found border at "+i);
		return i;
		}
		
		
	protected static int findHorCropPos(ImagePlus im, int start,int dir, double baseline)
		{
		double med;
		int w=im.getWidth();
		int h=im.getHeight();
		int i;
		int[] medValues;
		int right=w/2 - 10;
		int[] rights;
		double typical=PixelMean.get(im,right,0,20,h);
		double thresh=(typical+baseline)/2;
		medValues=new int[Math.round(w/4)+1];
		rights=new int[Math.round(w/4)+1];

		for (i=0;i<h/4;i++)
			{
			right=start+i*1*dir+(dir-1);
			rights[i]=right;
			med=PixelMean.get(im,right,0,2,h);
			if (med>thresh)
				{ medValues[i]=1; }
			else
				{ medValues[i]=0; }
			}
			
		i=i-2;
		int flag=1;
		while ((flag>0) & (i>8))
			{
			if (partialSum(medValues,i-9,i) == 0)
				{ flag=0; }	
			i=i-1;
			}
		i=rights[i];
		return i;
		}
	
	
	private void setCroppingParameters(int gridType)
		{
		switch (gridType)
			{
			case ColonyGrid.GRID_96:
			PLATE_EDGE_HEIGHT_FRAC=ColonyGrid96.PLATE_EDGE_HEIGHT_FRAC;
			PLATE_EDGE_WIDTH_FRAC=ColonyGrid96.PLATE_EDGE_WIDTH_FRAC;
			MINIMUM_PLATE_HEIGHT=ColonyGrid96.MINIMUM_PLATE_HEIGHT;
			MINIMUM_PLATE_WIDTH=ColonyGrid96.MINIMUM_PLATE_WIDTH;
			break;
			case ColonyGrid.GRID_384:
			PLATE_EDGE_HEIGHT_FRAC=ColonyGrid384.PLATE_EDGE_HEIGHT_FRAC;
			PLATE_EDGE_WIDTH_FRAC=ColonyGrid384.PLATE_EDGE_WIDTH_FRAC;
			MINIMUM_PLATE_HEIGHT=ColonyGrid384.MINIMUM_PLATE_HEIGHT;
			MINIMUM_PLATE_WIDTH=ColonyGrid384.MINIMUM_PLATE_WIDTH;
			break;
			case ColonyGrid.GRID_768a:
			PLATE_EDGE_HEIGHT_FRAC=ColonyGrid768a.PLATE_EDGE_HEIGHT_FRAC;
			PLATE_EDGE_WIDTH_FRAC=ColonyGrid768a.PLATE_EDGE_WIDTH_FRAC;
			MINIMUM_PLATE_HEIGHT=ColonyGrid768a.MINIMUM_PLATE_HEIGHT;
			MINIMUM_PLATE_WIDTH=ColonyGrid768a.MINIMUM_PLATE_WIDTH;
			break;
			}
		}
	}
