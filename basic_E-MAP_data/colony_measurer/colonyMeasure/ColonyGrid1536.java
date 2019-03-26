/*
 * Created on Aug 29, 2006
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

import java.awt.Rectangle;

import ij.ImagePlus;
import ij.gui.NewImage;
import ij.gui.OvalRoi;
import ij.measure.ResultsTablePlus;
import ij.process.ByteBlitter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

public class ColonyGrid1536 extends ColonyGrid 
	{
	
	static int SEARCH_WIDTH = 60;

	static final double PLATE_EDGE_HEIGHT_FRAC=0.030;
	static final double PLATE_EDGE_WIDTH_FRAC=0.03;
	static final int MINIMUM_PLATE_HEIGHT=1150;//1200;
	static final int MINIMUM_PLATE_WIDTH=1650;//1700;
	static final int NUM_ROW=32;
	static final int NUM_COL=48;
	static final int NUM_DUP=1;
	int imageHeight=0;

	public ColonyGrid1536() {}
	
	public ColonyGrid1536(ImagePlus im, ResultsTablePlus rt)
		{
		numRow = NUM_ROW;
		numCol = NUM_COL;
		numDup = NUM_DUP;
		heightTol = 200;
		super.setUpGrid(im, rt);
		}

		
	public ColonyGrid1536(ImagePlus im, int x1, int y1, int x2, int y2)
		{
		numRow = NUM_ROW;
		numCol = NUM_COL;
		numDup = NUM_DUP;
		heightTol = 200;
		setUpGrid(im, x1, y1, x2, y2);
		}


	public ColonyGrid1536(ImagePlus im, ResultsTablePlus rt, ColonyGrid old_grid, int rep)
		{
		//For taking a second shot at an image
		numRow = NUM_ROW;
		numCol = NUM_COL;
		numDup = NUM_DUP;
		heightTol = 200;
		super.setUpGrid(im, rt, old_grid, rep);
		}
		
		
	public void setUpGrid(ImagePlus im, int x1, int y1, int x2, int y2)
		{
		//This version is used for the manual processing
		slope = ((double)y2-(double)y1)/((double)x2-(double)x1);
		gridWidth = ((double)x2-(double)x1);
		fullSpace = gridWidth/47;
		gridHeight = fullSpace*31;
		leftBound = (double)x1+slope*gridHeight/2.0;
		gridTop = (double)y1+slope*gridWidth/2.0;
		halfSpace = fullSpace/2.0;
		gridBoxW = halfSpace*1.9;
		gridBoxD = gridBoxW/2;
		}
		
		
	protected void setGridParameters(ImagePlus im, ResultsTablePlus rt, int w, int h)
		{
		gridWidth = (rightBound-leftBound);
		fullSpace = gridWidth/47;
		halfSpace = fullSpace/2;
		gridHeight = fullSpace*31;
		gridTop = (h-gridHeight)/2;	//Initial guess
		gridBoxW = halfSpace*1.9;//
		gridBoxD=gridBoxW/2;
		}

	
	public OvalRoi getBoundingOval(int i, int j, int k)
		{
		//returns an oval that should enclose the colony at the grid position (i,j) (colony replicate k)
		int top;
		int left;
		int width = (int)Math.round(gridBoxW);
		int height = (int)Math.round(gridBoxW);
		
		//There is no second colony in this grid
		top = (int)Math.round(gridTop+j*fullSpace - (23.5-i)/47*gridWidth*slope - gridBoxD);
		left = (int)Math.round(leftBound-gridBoxD + i*fullSpace + (15.5-j)/31*gridHeight*slope);
		OvalRoi roi = new OvalRoi(left,top,width,height);
		
		return roi;
		}
	
	
	public OvalRoi getBoundingOvalForScoring(int i, int j, int k)
		{
		//returns an oval that should enclose the colony at the grid position (i,j) (colony replicate k)
		int top;
		int left;
		int width=(int)Math.round(gridBoxW);
		int height=(int)Math.round(gridBoxW);
		
		//There is no second colony in this grid
		top=(int)Math.round(gridTop+j*fullSpace - (23.5-i)/47*gridWidth*slope - gridBoxD);
		left=(int)Math.round(leftBound-gridBoxD + i*fullSpace + (15.5-j)/31*gridHeight*slope);
		
		OvalRoi roi=new OvalRoi(left,top + k*imageHeight,width,height);
		if (left < 0)
			{ System.out.println("Neg val"); }

		return roi;
		}	
	}
