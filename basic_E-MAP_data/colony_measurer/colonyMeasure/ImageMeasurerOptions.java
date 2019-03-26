/*
 * Created on Apr 3, 2006
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
public class ImageMeasurerOptions {
	int gridType;
	int cropType;
	int showOption;
	int preprocType;
	int threshold;
	
	// added by Assen
	int writeMaskedOption;
	
	public ImageMeasurerOptions(int gT, int cT, int sO, int wmO, int thres) {
		gridType=gT;
		cropType=cT;
		showOption=sO;
		writeMaskedOption = wmO;
		threshold = thres;
	}
}
