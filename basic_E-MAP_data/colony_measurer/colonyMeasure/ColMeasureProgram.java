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

import ij.ImagePlus;
import ij.process.*;
import ij.gui.ImageCanvas;
import ij.gui.Line;
import ij.gui.PointRoi;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

// Added by Assen
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.*;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;

public class ColMeasureProgram extends JFrame implements ActionListener, ItemListener {

	JMenuItem itemMeasure, itemMeasureFolder, itemManualMeasure, itemQuit, itemTimeStackMeasure;
	JMenuItem itemShowAbout;
	JCheckBoxMenuItem itemGrid96, itemGrid384, itemGrid768a, itemGrid768b, itemGrid1536;
	JCheckBoxMenuItem itemAutoCropDecide, itemAutoCrop, itemPreCrop;

	// added by Assen
	JCheckBoxMenuItem itemSaveGridOverlayImgFile, itemSaveGridOverlayImgFolder, itemSaveGridOverlayImgManual;
	
	
	JMenuBar mb;
	JPanel cp;
	JFileChooser fc;
	JLabel statusLabel;
	int numProcessed;
	int numErrorsFound;
	int gridType;
	int cropType;
	String errorFileNames[];
	int errorCodes[];
	
	static int ERR_NONE = 0;
	static int ERR_OTHER = 1;
	static int ERR_TOO_FEW_COL = 2;
	static int ERR_TOO_MANY_COL = 3;
	static int ERR_CROP_ERROR = 4;
	static int ERR_REPLICATES = 5;
	
	// 0 for threshold calc only, 1 for full analysis
	static int AN_THRES_ONLY = 0;
	static int AN_FULL = 1;
	
	static int ZERO_THRES = 0;		// assen: silly but makes me understand it better
	
	// defaults for saving the image-grid overlay, added by Assen
	static int WRITE_MASKED_OPTION_FILE = 0;
	static int WRITE_MASKED_OPTION_FOLDER = 1;
	static int WRITE_MASKED_OPTION_MANUAL = 0;
	
	// error messages
	static String[] ERR_MESSAGE = {
		"",
		"Processing required multiple attempts - may be incorrect:",
		"Grid Alignment - Too few colonies:",
		"Grid Alignment - Too many colonies:",
		"Possible cropping error:",
		"Replicate colony measurements are inconsistent:"
		}; 
	
	
	public static void main(String[] args)
		{ ColMeasureProgram cmp = new ColMeasureProgram(); }
	
	public ColMeasureProgram()
		{
		numErrorsFound=0;

		//Default values
		gridType = ColonyGrid.GRID_384;
		cropType = ColonyGrid.AUTO_CROP_OPTION;
				
		setTitle("Colony Measuring Program");
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		cp = new JPanel(new BorderLayout());
		setContentPane(cp);

		// added by assen
		cp.setSize(400, 50);
		
		// menubar stuff
		mb = new JMenuBar();
		mb.setBackground(Color.cyan);
		buildFileMenu();
		buildGridTypeMenu();
		buildCropTypeMenu();
		
		// added by assen
		buildOptionsMenu();
		buildAboutMenu();
		
		statusLabel=new JLabel("Ready");
		cp.add(statusLabel,BorderLayout.NORTH);
		initFileChooser();
		setJMenuBar(mb);
		mb.validate();
		setVisible(true);
		pack();

		// added by Assen
		setSize(400,100);
		setResizable(false);

		validate();		
		}


	private void buildGridTypeMenu()
		{
		JMenu gridMen;
		
		gridMen = new JMenu("Grid    ");
		gridMen.setBackground(Color.blue);
		itemGrid96 = new JCheckBoxMenuItem("96 colonies",false);
		itemGrid96.addActionListener(this);
		gridMen.add(itemGrid96);
		
		// Assen
		itemGrid384 = new JCheckBoxMenuItem("384 colonies",true);
		itemGrid384.addActionListener(this);
		gridMen.add(itemGrid384);
		

		itemGrid768a = new JCheckBoxMenuItem("768 colonies - diagonal replicates",false);
		itemGrid768a.addActionListener(this);
		gridMen.add(itemGrid768a);
		
		itemGrid1536 = new JCheckBoxMenuItem("1536 colonies",false);
		itemGrid1536.addActionListener(this);
		gridMen.add(itemGrid1536);


		mb.add(gridMen);
		}
		

	// added by Assen	
	private void buildOptionsMenu()
		{
		JMenu optionsMen;
			
		optionsMen = new JMenu("Options    ");
		optionsMen.setBackground(Color.green);
		
		// items
		itemSaveGridOverlayImgFile = new JCheckBoxMenuItem("Save Grid Image Single", false);
		itemSaveGridOverlayImgFile.addItemListener(this);
		optionsMen.add(itemSaveGridOverlayImgFile);

		itemSaveGridOverlayImgFolder = new JCheckBoxMenuItem("Save Grid Image Folder and Timestack", true);
		itemSaveGridOverlayImgFolder.addItemListener(this);
		optionsMen.add(itemSaveGridOverlayImgFolder);


		itemSaveGridOverlayImgManual = new JCheckBoxMenuItem("Save Grid Image Manual", false);
		itemSaveGridOverlayImgManual.addItemListener(this);
		optionsMen.add(itemSaveGridOverlayImgManual);
		
		// add to the main menu
		mb.add(optionsMen);
		}	
		
	private void buildCropTypeMenu()
		{
		JMenu cropMen;

		cropMen = new JMenu("Cropping    ");
		cropMen.setBackground(Color.red);

		itemAutoCropDecide = new JCheckBoxMenuItem("Automatically Choose Method",false);
		itemAutoCropDecide.addActionListener(this);
		cropMen.add(itemAutoCropDecide);

		itemAutoCrop = new JCheckBoxMenuItem("Always Autodetect Plate Edges",true);
		itemAutoCrop.addActionListener(this);
		cropMen.add(itemAutoCrop);

		itemPreCrop = new JCheckBoxMenuItem("Images Are Already Cropped to Plate Edges",false);
		itemPreCrop.addActionListener(this);
		cropMen.add(itemPreCrop);

		mb.add(cropMen);
		}	
	
	private void buildFileMenu() {
		JMenu fileMen;
		
		fileMen = new JMenu("File    ");
		fileMen.setBackground(Color.cyan);

		itemMeasure = new JMenuItem("Single file mode");
		itemMeasure.addActionListener(this);
		fileMen.add(itemMeasure);

		itemMeasureFolder = new JMenuItem("Batch mode");
		itemMeasureFolder.addActionListener(this);
		fileMen.add(itemMeasureFolder);

		itemManualMeasure = new JMenuItem("Manual mode");
		itemManualMeasure.addActionListener(this);
		fileMen.add(itemManualMeasure);
		
		// added by assen
		itemTimeStackMeasure = new JMenuItem("Time stack mode");
		itemTimeStackMeasure.addActionListener(this);
		fileMen.add(itemTimeStackMeasure);


		itemQuit = new JMenuItem("Exit");
		itemQuit.addActionListener(this);
		fileMen.add(itemQuit);

		mb.add(fileMen);
	}
	
	private void buildAboutMenu() {
		JMenu aboutMen;

		aboutMen = new JMenu("About");
		aboutMen.setBackground(Color.cyan);

		itemShowAbout = new JMenuItem("Show info");
		itemShowAbout.addActionListener(this);
		aboutMen.add(itemShowAbout);

		mb.add(aboutMen);
	}

	
	private void initFileChooser()
		{
		fc = new JFileChooser();
		//String sepchar = File.separator;
		//String rootDirName = "C:"+File.separator;
		//File rootDir = new File(rootDirName);
		//fc.setCurrentDirectory(rootDir);
		}
	
	
	// assen changed header from void to int
	private int processFile(String path, int showOption, int wmO, int m, int threshold)
		{
		int t = 0;
		String temp[];
		int temp2[];
		int i;
		int mode = m;

		if ((path.endsWith(".jpg")) | (path.endsWith(".JPG")))
			{
			System.out.println("ColMeasureProgram->processFile input: "+path);	
			numProcessed++;
			ImagePlus im = new ImagePlus(path);
			ImageProcessor ip = im.getProcessor();
			String filename=path+".dat";
			System.out.println("ColMeasureProgram->processFile output: "+filename); 

			ImageMeasurerOptions imo = new ImageMeasurerOptions(gridType, cropType, showOption, wmO, threshold);
			ImageMeasurer imM = new ImageMeasurer(ip, filename, imo);

			// check to see if we have a precalculated threshold and set it in the imM object
			if ( (mode == AN_FULL) && (threshold != ZERO_THRES) ) {imM.THRESHOLD = threshold; }
			imM.fullAnalysis(mode);
			
			t = imM.THRESHOLD;
						
			if (showOption==ImageMeasurer.SHOW_NOTHING)
				{ updateLabel(statusLabel,"Processed " + numProcessed + " files"); }
				
			if (imM.ERROR_FLAG>0)
				{
				temp=new String[numErrorsFound+1];
				temp2=new int[numErrorsFound+1];
				
				for (i=0;i<numErrorsFound;i++)
					{
					temp[i]=errorFileNames[i];
					temp2[i]=errorCodes[i];
					}
					
				temp[numErrorsFound]=path;
				temp2[numErrorsFound]=imM.ERROR_FLAG;
				numErrorsFound++;
				errorFileNames=temp;
				errorCodes=temp2;
				}				
			}
		return t;		
		}

// added by Assen
	private void processTimeStack(File file_time0)
		{	
		int thres0;	
		// process time0
		System.out.println("ColMeasureProgram->processTimeStack starting image: "+file_time0.getPath()); 	
		thres0 = processFile(file_time0.getPath(), ImageMeasurer.SHOW_NOTHING, WRITE_MASKED_OPTION_FOLDER, AN_THRES_ONLY, ZERO_THRES);
		
		// extract directory and process all files there doing full analysis
		File parent_dir = file_time0.getParentFile();
		
		// System.out.println("processTimeStack Dir: " + parent_dir);
		processFolder(parent_dir, thres0);
		}
		
		
	private void processFolder(File dir, int thres)
		{
		System.out.println("Batch mode Dir: " + dir);	
		File[] children = dir.listFiles(); 
		if (children == null)
			{ int t = processFile(dir.getPath(), ImageMeasurer.SHOW_NOTHING, WRITE_MASKED_OPTION_FOLDER, AN_FULL, thres); } 
		else
			{
			for (int i=0; i < children.length; i++)
				{ processFolder(children[i], thres); }
			}
		}
		
		
	private void processFileManual(String path, int showOption)
		{
		if ((path.endsWith(".jpg"))|(path.endsWith(".JPG")))
			{
			numProcessed=numProcessed+1;
			ImagePlus im = new ImagePlus(path);
			ImageProcessor ip = im.getProcessor();
			String filename=path+".dat";
			ManualProcessor mp=new ManualProcessor(ip,filename,this);
			}
		}
		
//	added by Assen	
	public void itemStateChanged(ItemEvent e)
		{	
		if (e.getSource() == itemSaveGridOverlayImgFile)
			{
			if (itemSaveGridOverlayImgFile.getState() == true)
				{	WRITE_MASKED_OPTION_FILE = 1; }
			else
				{	WRITE_MASKED_OPTION_FILE = 0; }
			}
			
		if (e.getSource() == itemSaveGridOverlayImgFolder)
			{
			if (itemSaveGridOverlayImgFolder.getState() == true)
				{	WRITE_MASKED_OPTION_FOLDER = 1; }
			else 
				{	WRITE_MASKED_OPTION_FOLDER = 0;	}
			}
		
		if (e.getSource() == itemSaveGridOverlayImgManual)
			{
			if (itemSaveGridOverlayImgManual.getState() == true)
				{	WRITE_MASKED_OPTION_MANUAL = 0; }
			else
				{	WRITE_MASKED_OPTION_MANUAL = 1;	}
			}
		}		


// action listener, things will get messy here	
	public void actionPerformed(ActionEvent e)
		{
		File selectedFile;
	
	// single file mode	
		if (e.getSource() == itemMeasure)
			{
			fc.setFileSelectionMode(fc.FILES_ONLY);
			int returnval = fc.showOpenDialog(this);
			
			if (returnval == JFileChooser.APPROVE_OPTION)
				{
				selectedFile = fc.getSelectedFile();
				String path = selectedFile.getPath();
				numErrorsFound=0;
				 
				int t = processFile(path,ImageMeasurer.SHOW_GRID_OVERLAY, WRITE_MASKED_OPTION_FILE, AN_FULL, ZERO_THRES);
				updateLabel(statusLabel,"Done processing image");
				
				if (numErrorsFound > 0)
					{
					updateLabel(statusLabel,"There may have been an error processing this file.");
					File dir = selectedFile.getParentFile();
					exportErrorList(dir);
					}
				}
			}

	// manual mode
		if (e.getSource() == itemManualMeasure)
			{
			fc.setFileSelectionMode(fc.FILES_ONLY);
			int returnval = fc.showOpenDialog(this);

			if (returnval == JFileChooser.APPROVE_OPTION)
				{
				selectedFile = fc.getSelectedFile();
				String path = selectedFile.getPath();
				numErrorsFound=0;
				processFileManual(path,ImageMeasurer.SHOW_GRID_OVERLAY);
				updateLabel(statusLabel,"Done processing image");
	
				if (numErrorsFound>0)
					{
					updateLabel(statusLabel,"There may have been an error processing this file.");
					File dir=selectedFile.getParentFile();
					exportErrorList(dir);
					}
				}
			}
	
	// batch mode		
		if (e.getSource() == itemMeasureFolder)
			{
			fc.setFileSelectionMode(fc.DIRECTORIES_ONLY);
			int returnval = fc.showOpenDialog(this);
			
			if (returnval == JFileChooser.APPROVE_OPTION)
				{
				File dir = fc.getSelectedFile();
				numProcessed=0;
				numErrorsFound=0;
				processFolder(dir, ZERO_THRES);
				updateLabel(statusLabel,"Done processing folder");
				
				if (numErrorsFound>0)
					{ updateLabel(statusLabel,"There may have been an error processing at least one of the files. Dumping report."); exportErrorList(dir); }
				}
			}
			

	// timestack mode
	 	if (e.getSource() == itemTimeStackMeasure)
			{
			fc.setFileSelectionMode(fc.FILES_ONLY);
			int returnval = fc.showOpenDialog(this);
			
			if (returnval == JFileChooser.APPROVE_OPTION)
				{
				selectedFile = fc.getSelectedFile();
				numErrorsFound=0;
				 
				processTimeStack(selectedFile);
				updateLabel(statusLabel,"Done processing timestack");
				
				if (numErrorsFound>0)
					{ 
					updateLabel(statusLabel,"There may have been an error processing this file."); 
					// exportErrorList(dir);
					}
				}
			}	
	
			
	// exit, dunno if tis needed, the setDefaultCloseOperation(EXIT_ON_CLOSE) up there should do the job		
		if (e.getSource() == itemQuit)
			{	this.dispose();	}
			
	
	// grid 96		
		if (e.getSource() == itemGrid96)
			{
			gridType=ColonyGrid.GRID_96;
			itemGrid96.setState(true);
			itemGrid384.setState(false);
			itemGrid768a.setState(false);
			itemGrid1536.setState(false);
			//itemGrid768b.setState(false);
			}

	// grid 384
		if (e.getSource() == itemGrid384)
			{
			gridType=ColonyGrid.GRID_384;
			itemGrid96.setState(false);
			itemGrid384.setState(true);
			itemGrid768a.setState(false);
			itemGrid1536.setState(false);
			//itemGrid768b.setState(false);
			}

	// grid 768
		if (e.getSource() == itemGrid768a)
			{
			gridType=ColonyGrid.GRID_768a;
			itemGrid96.setState(false);
			itemGrid384.setState(false);
			itemGrid768a.setState(true);
			itemGrid1536.setState(false);
			//itemGrid768b.setState(false);
			}
	
	// grid 1536		
		if (e.getSource() == itemGrid1536)
			{
			gridType=ColonyGrid.GRID_1536;
			itemGrid96.setState(false);
			itemGrid384.setState(false);
			itemGrid768a.setState(false);
			itemGrid1536.setState(true);
			//itemGrid768b.setState(false);
			}
			

	// cropping stuff
		if (e.getSource() == itemAutoCropDecide)
			{
			cropType=ColonyGrid.AUTO_CROP_DECIDE_OPTION;
			itemAutoCropDecide.setState(true);
			itemAutoCrop.setState(false);
			itemPreCrop.setState(false);
			}

		if (e.getSource() == itemAutoCrop)
			{
			cropType=ColonyGrid.AUTO_CROP_OPTION;
			itemAutoCropDecide.setState(false);
			itemAutoCrop.setState(true);
			itemPreCrop.setState(false);
			}

		if (e.getSource() == itemPreCrop)
			{
			cropType=ColonyGrid.PRECROP_OPTION;
			itemAutoCropDecide.setState(false);
			itemAutoCrop.setState(false);
			itemPreCrop.setState(true);
			}

	// the about menu
		if (e.getSource() == itemShowAbout)
			{
			ColMeasureAboutWindow cma=new ColMeasureAboutWindow();
			}
		}
		
	private void writeln(BufferedWriter bw, String str)
		{
		try
			{
			bw.write(str, 0, str.length());
			bw.newLine();
			}
			
		catch (IOException io)
			{	System.out.println("Error writing a line to the file");	}
		}
	
	
	public void exportErrorList(File dir)
		{
		BufferedWriter bw;
		File outFile;
		int i,j;
		String filename="errors.txt";

		outFile=new File(dir,filename);
		try
			{
			bw = new BufferedWriter(new FileWriter(outFile));
			writeln(bw,"Files for which the processing may not have worked correctly");
			for (i=5;i>0;i--)
				{ //iterate over the error types
				writeln(bw,"");
				writeln(bw,ERR_MESSAGE[i]);
				
				for (j=0;j<errorFileNames.length;j++)
					{
					if (errorCodes[j]==i)
						{	writeln(bw,errorFileNames[j]);	}
					}
				}
			bw.close();
			}
			
		catch (FileNotFoundException fnf)
			{	System.err.println("ColMeasureProgram->exportErrorList: Could not make the error file");	}

		catch (IOException io)
			{	System.err.println("ColMeasureProgram->exportErrorList: Error reading the file");	}		
		}


	private void updateLabel(JLabel label, String str)
		{
		label.setText(str);
		update(this.getGraphics());
		}
	}
//#################### the end of it
