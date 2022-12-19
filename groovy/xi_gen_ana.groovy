import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.clas.physics.RecEvent;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.ui.TCanvas;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.DataLine;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.*;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.detector.DetectorResponse;
import org.jlab.clas.detector.DetectorParticle;
import org.jlab.detector.base.DetectorType;
import org.jlab.service.ec.*;
import org.jlab.geom.prim.Vector3D;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import javax.swing.JFrame;
import java.util.Random;
//import org.apache.commons.math3.distribution;


GStyle.setPalette("kRainBow");
GStyle.getH1FAttributes().setOptStat("1110");

//double ebeam = 10.604;
double ebeam = 10.6;
final int BLUE = 9;
final int LIGHTGREEN = 3;
	
	H2F hi_q2_w = new H2F("hi_q2_w", "hi_q2_w", 100, 3, 5, 100, 0.00001, 0.5);
	hi_q2_w.setTitleX("W (GeV)");
	hi_q2_w.setTitleY("Q^2 (GeV^2)")

	H1F hi_w = new H1F("hi_w", "hi_w", 100, 3, 5);
	hi_w.setTitleX("W (GeV)");
	hi_w.setFillColor(LIGHTGREEN);

	H1F hi_q2 = new H1F("hi_q2", "hi_q2", 100, 0.00001, 0.5);
	hi_q2.setTitleX("Q^2 (GeV^2)");
	hi_q2.setFillColor(LIGHTGREEN);

	// generated eletron
	H2F hi_FT_e_p_the = new H2F("hi_FT_e_p_the", "hi_FT_e_p_the", 100, 0, 6, 100, 0, 8);		
	hi_FT_e_p_the.setTitle("electron p vs #theta (^o)");
	hi_FT_e_p_the.setTitleX("#theta (^o)");
	hi_FT_e_p_the.setTitleY("p (GeV)");

	H2F hi_FT_e_p_f = new H2F("hi_FT_e_p_f", "hi_FT_e_p_f", 100, -180, 180, 100, 0, 8);
	hi_FT_e_p_f.setTitle("electron p vs #phi");
	hi_FT_e_p_f.setTitleX("#phi (^o)");
	hi_FT_e_p_f.setTitleY("p (GeV)");

	H2F hi_FT_e_the_phi = new H2F("hi_FT_e_the_phi", "hi_FT_e_the_phi", 100, -180, 180, 100, 0, 6);
	hi_FT_e_the_phi.setTitle("electron #theta vs #phi");
	hi_FT_e_the_phi.setTitleX("#phi (^o)");
	hi_FT_e_the_phi.setTitleY("#theta (^o)");


	H1F hi_e_vz = new H1F("hi_e_vz", "hi_e_vz", 200, -20, 20);
	hi_e_vz.setTitleX("vz");
	hi_e_vz.setFillColor(LIGHTGREEN);

	H1F hi_e_vx = new H1F("hi_e_vx", "hi_e_vx", 200, -20, 20);
	hi_e_vx.setTitleX("vx");
	hi_e_vx.setFillColor(LIGHTGREEN);

	H1F hi_e_vy = new H1F("hi_e_vy", "hi_e_vy", 200, -20, 20);
	hi_e_vy.setTitleX("vy");
	hi_e_vy.setFillColor(LIGHTGREEN);

	H1F hi_virphoton = new H1F("hi_virphoton", "hi_virphoton", 100, 0, 12);
	hi_virphoton.setTitleX("E_#gamma (GeV)");
	hi_virphoton.setFillColor(LIGHTGREEN);

	// generated kps

	// KP1
	H2F hi_kp1_p_the = new H2F("hi_kp1_p_the", "hi_kp1_p_the", 100, 0, 100, 100, 0, 8);
	hi_kp1_p_the.setTitleX("#theta");
	hi_kp1_p_the.setTitleY("p (GeV)");

	H2F hi_kp1_p_phi = new H2F("hi_kp1_p_phi", "hi_kp1_p_phi", 100, -180, 180, 100, 0, 8);
	hi_kp1_p_phi.setTitleX("#phi (^o)");
	hi_kp1_p_phi.setTitleY("p (GeV)");

	H1F hi_kp1_vz = new H1F("hi_kp1_vz", "hi_kp1_vz", 200, -20, 20);
	hi_kp1_vz.setFillColor(LIGHTGREEN);

	// KP2
	H2F hi_kp2_p_the = new H2F("hi_kp2_p_the", "hi_kp2_p_the", 100, 0, 100, 100, 0, 8);
	hi_kp2_p_the.setTitleX("#theta (^o)");
	hi_kp2_p_the.setTitleY("p (GeV)");

	H2F hi_kp2_p_phi = new H2F("hi_kp2_p_phi", "hi_kp2_p_phi", 100, -180, 100, 100, 0, 8);
	hi_kp2_p_phi.setTitleX("#phi (^o)");
	hi_kp2_p_phi.setTitleY("p (GeV)");

	H1F hi_kp2_vz = new H1F("hi_kp2_vz", "hi_kp2_vz", 200, -20, 20);
	hi_kp2_vz.setFillColor(LIGHTGREEN);
	//KM
	H2F hi_km_p_the = new H2F("hi_km_p_the", "hi_km_p_the", 100, 0, 100, 100, 0, 8);
	hi_km_p_the.setTitleX("#theta");
	hi_km_p_the.setTitleY("p (GeV)");

	H2F hi_km_p_phi = new H2F("hi_km_p_phi", "hi_km_p_phi", 100, -180, 180, 100, 0, 8);
	hi_km_p_phi.setTitleX("#phi (^o)");
	hi_km_p_phi.setTitleY("p (GeV)");

	H1F hi_km_vz = new H1F("hi_km_vz", "hi_km_vz", 200, -20, 20);
	hi_km_vz.setFillColor(LIGHTGREEN);
	//PROton
	H2F hi_pro_p_the = new H2F("hi_pro_p_the", "hi_pro_p_the", 100, 0, 100, 100, 0, 8);
	hi_pro_p_the.setTitleX("#theta (^o)");
	hi_pro_p_the.setTitleY("p (GeV)");

	H2F hi_pro_p_phi = new H2F("hi_pro_p_phi", "hi_pro_p_phi", 100, -180, 180, 100, 0, 8);
	hi_pro_p_phi.setTitleX("#phi (^o)");
	hi_pro_p_phi.setTitleY("p (GeV)");

	H1F hi_pro_vz = new H1F("hi_pro_vz", "hi_pro_vz", 200, -20, 20);
	hi_pro_vz.setFillColor(LIGHTGREEN);
	//PIM
	H2F hi_pim_p_the = new H2F("hi_pim_p_the", "hi_pim_p_the", 100, 0, 100, 100, 0, 8);
	hi_pim_p_the.setTitleX("#theta (^o)");
	hi_pim_p_the.setTitleY("p (GeV)");

	H2F hi_pim_p_phi = new H2F("hi_pim_p_phi", "hi_pim_p_phi", 100, -180, 180, 100, 0, 8);
	hi_pim_p_phi.setTitleX("#phi (^o)");
	hi_pim_p_phi.setTitleY("p (GeV)");

	H1F hi_pim_vz = new H1F("hi_pim_vz", "hi_pim_vz", 200, -20, 20);
	hi_pim_vz.setFillColor(LIGHTGREEN);

	H1F hi_mm_ekpkpkm = new H1F("hi_mm_ekpkpkm", "hi_mm_ekpkpkm", 50, 0.9, 1.4);
	hi_mm_ekpkpkm.setFillColor(LIGHTGREEN);

	F1D f1_xi = new F1D("f1_xi", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
	f1_xi.setParameter(0, 0);
    	f1_xi.setParameter(1, 1);
    	f1_xi.setParameter(2, 0.2);
    	f1_xi.setLineWidth(2);
    	f1_xi.setLineColor(2);
    	f1_xi.setOptStat("1111");

	H1F hi_mm_ekpkp   = new H1F("hi_mm_ekpkp", "hi_mm_ekpkp", 50, 1.6, 3.1);//25, 0.0, 2.4//50, 1.6, 2.1
	hi_mm_ekpkp.setFillColor(LIGHTGREEN);

	DataGroup dg_mm = new DataGroup(2,1);
	dg_mm.addDataSet(hi_mm_ekpkp, 0);
	dg_mm.addDataSet(f1_xi, 0);
	dg_mm.addDataSet(hi_mm_ekpkpkm, 1);


	DataGroup dg_electron = new DataGroup(3,4);
	dg_electron.addDataSet(hi_FT_e_p_the, 0);
	dg_electron.addDataSet(hi_FT_e_p_f, 1);
	dg_electron.addDataSet(hi_FT_e_the_phi, 2);
	dg_electron.addDataSet(hi_e_vx, 3);
	dg_electron.addDataSet(hi_e_vy, 4);
	dg_electron.addDataSet(hi_e_vz, 5);
	dg_electron.addDataSet(hi_q2, 6);
	dg_electron.addDataSet(hi_w, 7);
	dg_electron.addDataSet(hi_virphoton, 8);
	dg_electron.addDataSet(hi_q2_w, 9);

	DataGroup dg_kp1 = new DataGroup(2, 2);
	dg_kp1.addDataSet(hi_kp1_p_the, 0);
	dg_kp1.addDataSet(hi_kp1_p_phi, 1);
	dg_kp1.addDataSet(hi_kp1_vz, 2);

	DataGroup dg_kp2 = new DataGroup(2, 2);
	dg_kp2.addDataSet(hi_kp2_p_the, 0);
	dg_kp2.addDataSet(hi_kp2_p_phi, 1);
	dg_kp2.addDataSet(hi_kp2_vz, 2);

	DataGroup dg_km = new DataGroup(2,2);
	dg_km.addDataSet(hi_km_p_the, 0);
	dg_km.addDataSet(hi_km_p_phi, 1);
	dg_km.addDataSet(hi_km_vz, 2);


	DataGroup dg_pro = new DataGroup(2,2);
	dg_pro.addDataSet(hi_pro_p_the, 0);
	dg_pro.addDataSet(hi_pro_p_phi, 1);
	dg_pro.addDataSet(hi_pro_vz, 2);

	DataGroup dg_pim = new DataGroup(2,2);
	dg_pim.addDataSet(hi_pim_p_the, 0);
	dg_pim.addDataSet(hi_pim_p_phi, 1);
	dg_pim.addDataSet(hi_pim_vz, 2);


int nevent = -1;
for(int k = 0; k < args.length; k++) {


	HipoDataSource reader = new HipoDataSource();
    reader.open(args[k]);

    while(reader.hasEvent() == true && nevent < 10000000) {
    	DataEvent event = reader.getNextEvent();
    	nevent++;
    	if (nevent%100000 == 0) System.out.println("Analyzed " + nevent + " events");
    	//event.show();

    	DataBank mcParticle = null;

    	if(event.hasBank("MC::Particle"))			mcParticle =  event.getBank("MC::Particle");

    	Particle mcEl = null;
    	Particle mcKp1 = null;
    	Particle mcKp2 = null;
    	Particle mcKm  = null;
    	Particle mcPro = null;
    	Particle mcPim = null;

    	LorentzVector virtualPhoton  = null;
   		LorentzVector hadronSystem   = null;
   		LorentzVector virtualPhotonP = null;
   		LorentzVector hadronSystemP  = null;
   		LorentzVector xi		 = null;
   		LorentzVector lambda          = null;

    	if(event.hasBank("MC::Particle") == true ){
    		
    		List<Particle> mckps = new ArrayList<Particle>();

    		DataBank mcBank = event.getBank("MC::Particle");

    		for(int loop = 0; loop < mcBank.rows(); loop++) {

    			//NormalDistribution smFactor = new NormalDistribution(0, 0.01);
    			Random rand = new Random();
    			double smearFactor = rand.nextGaussian()/100;
    			System.out.println("Momentum Smearing factor for MC" + " " + smearFactor);

				if(mcBank.getInt("pid", loop) == 11) {
					mcEl = new Particle( 
				 					 mcBank.getInt("pid", loop), 
				 					 mcBank.getFloat("px", loop), 
				 					 mcBank.getFloat("py", loop), 
				 					 mcBank.getFloat("pz", loop),
				 					 mcBank.getFloat("vx", loop),
				 					 mcBank.getFloat("vy", loop),
				 					 mcBank.getFloat("vz", loop)); 
				}

				else if(mcBank.getInt("pid", loop) == 321) {

					Particle mcparticle = new Particle( 
				 					 mcBank.getInt("pid", loop), 
				 					 mcBank.getFloat("px", loop), 
				 					 mcBank.getFloat("py", loop), 
				 					 mcBank.getFloat("pz", loop),
				 					 mcBank.getFloat("vx", loop),
				 					 mcBank.getFloat("vy", loop),
				 					 mcBank.getFloat("vz", loop)); 

					double smearedP = mcparticle.p()*(1+smearFactor);
					mcparticle.setP(smearedP);
					mckps.add(mcparticle);
				}

				else if(mcBank.getInt("pid", loop) == -321) {
					mcKm = new Particle( 
				 					 mcBank.getInt("pid", loop), 
				 					 mcBank.getFloat("px", loop), 
				 					 mcBank.getFloat("py", loop), 
				 					 mcBank.getFloat("pz", loop),
				 					 mcBank.getFloat("vx", loop),
				 					 mcBank.getFloat("vy", loop),
				 					 mcBank.getFloat("vz", loop)); 
					double smearedP = mcKm.p()*(1+smearFactor);
					mcKm.setP(smearedP);
					//mckps.add(mcparticle);
				}

				else if(mcBank.getInt("pid", loop) == 2212) {
					mcPro = new Particle( 
				 					 mcBank.getInt("pid", loop), 
				 					 mcBank.getFloat("px", loop), 
				 					 mcBank.getFloat("py", loop), 
				 					 mcBank.getFloat("pz", loop),
				 					 mcBank.getFloat("vx", loop),
				 					 mcBank.getFloat("vy", loop),
				 					 mcBank.getFloat("vz", loop)); 
				}

				else if(mcBank.getInt("pid", loop) == -211) {
					mcPim = new Particle( 
				 					 mcBank.getInt("pid", loop), 
				 					 mcBank.getFloat("px", loop), 
				 					 mcBank.getFloat("py", loop), 
				 					 mcBank.getFloat("pz", loop),
				 					 mcBank.getFloat("vx", loop),
				 					 mcBank.getFloat("vy", loop),
				 					 mcBank.getFloat("vz", loop)); 
				}


    		}

    		if (mcEl != null) {

    			virtualPhoton = new LorentzVector(0.0, 0.0, ebeam, ebeam);
          		virtualPhoton.sub(mcEl.vector());
          		hadronSystem = new LorentzVector(0.0, 0.0, ebeam, 0.9383+ebeam);
          		hadronSystem.sub(mcEl.vector());
          		if(Math.toDegrees(mcEl.theta())>0){
          			dg_electron.getH1F("hi_q2").fill(-virtualPhoton.mass2());
          			dg_electron.getH1F("hi_w").fill(hadronSystem.mass());
          			dg_electron.getH1F("hi_e_vz").fill(mcEl.vz());
          			dg_electron.getH1F("hi_e_vx").fill(mcEl.vx());
          			dg_electron.getH1F("hi_e_vy").fill(mcEl.vy());
          			dg_electron.getH2F("hi_q2_w").fill(hadronSystem.mass(),-virtualPhoton.mass2());
          			dg_electron.getH2F("hi_FT_e_the_phi").fill(Math.toDegrees(mcEl.phi()), Math.toDegrees(mcEl.theta()));
    				dg_electron.getH2F("hi_FT_e_p_the").fill(Math.toDegrees(mcEl.theta()), mcEl.p());
    				dg_electron.getH2F("hi_FT_e_p_f").fill(Math.toDegrees(mcEl.phi()), mcEl.p());
    				dg_electron.getH1F("hi_virphoton").fill(ebeam - mcEl.e());
    			}
    		}

    		if (mckps.size() == 2 && mcEl != null && mcKm != null ){ //&& mcPro != null && mcPim != null
    			if(mckps.get(0).p() > mckps.get(1).p()){
    				mcKp1 = mckps.get(0);
    				mcKp2 = mckps.get(1);
    			} else {
    				mcKp1 = mckps.get(1);
    				mcKp2 = mckps.get(0);
    			}
    			
    			dg_kp1.getH2F("hi_kp1_p_the").fill(Math.toDegrees(mcKp1.theta()), mcKp1.p());
    			dg_kp1.getH2F("hi_kp1_p_phi").fill(Math.toDegrees(mcKp1.phi()), mcKp1.p());
    			dg_kp1.getH1F("hi_kp1_vz").fill(mcKp1.vz());
    			dg_kp2.getH2F("hi_kp2_p_the").fill(Math.toDegrees(mcKp2.theta()), mcKp2.p());
    			dg_kp2.getH2F("hi_kp2_p_phi").fill(Math.toDegrees(mcKp2.phi()), mcKp2.p());
    			dg_kp2.getH1F("hi_kp2_vz").fill(mcKp2.vz());


    			xi = new LorentzVector(0.0, 0.0, ebeam, 0.9383+ebeam);
    			xi.sub(mcEl.vector());
    			xi.sub(mcKp1.vector());
    			xi.sub(mcKp2.vector());

    			lambda = new LorentzVector(0.0, 0.0, ebeam, 0.9383+ebeam);
    			lambda.sub(mcEl.vector());
    			lambda.sub(mcKp1.vector());
    			lambda.sub(mcKp2.vector());
    			lambda.sub(mcKm.vector());

    			dg_mm.getH1F("hi_mm_ekpkp").fill(xi.mass());
    			dg_mm.getH1F("hi_mm_ekpkpkm").fill(lambda.mass()); 
    		}

    		if (mcKm != null) {

    			dg_km.getH2F("hi_km_p_the").fill(Math.toDegrees(mcKm.theta()), mcKm.p());
    			dg_km.getH2F("hi_km_p_phi").fill(Math.toDegrees(mcKm.phi()), mcKm.p());
    			dg_km.getH1F("hi_km_vz").fill(mcKm.vz());
    		}

    		if (mcPro != null) {
    			dg_pro.getH2F("hi_pro_p_the").fill(Math.toDegrees(mcPro.theta()), mcPro.p());
    			dg_pro.getH2F("hi_pro_p_phi").fill(Math.toDegrees(mcPro.phi()), mcPro.p());
    			dg_pro.getH1F("hi_pro_vz").fill(mcPro.vz());
    		}

    		if (mcPim != null) {
    			dg_pim.getH2F("hi_pim_p_the").fill(Math.toDegrees(mcPim.theta()), mcPim.p());
    			dg_pim.getH2F("hi_pim_p_phi").fill(Math.toDegrees(mcPim.phi()), mcPim.p());
    			dg_pim.getH1F("hi_pim_vz").fill(mcPim.vz());
    		} 
    	}
    }
}

//fitxi(dg_mm.getH1F("hi_mm_ekpkp"), dg_mm.getF1D("f1_xi"));

EmbeddedCanvasTabbed myCanvas = new EmbeddedCanvasTabbed("genElectron", "genKp1", "genKp2", "genKm", "genProton", "genPim", "MM");
myCanvas.getCanvas("genElectron").divide(3,4);
myCanvas.getCanvas("genElectron").setGridX(false);
myCanvas.getCanvas("genElectron").setGridY(false);
myCanvas.getCanvas("genElectron").setAxisFontSize(18);
myCanvas.getCanvas("genElectron").setAxisTitleSize(24);
myCanvas.getCanvas("genElectron").draw(dg_electron);
myCanvas.getCanvas("genElectron").getPad(0).getAxisZ().setLog(true);
myCanvas.getCanvas("genElectron").getPad(1).getAxisZ().setLog(true);
myCanvas.getCanvas("genElectron").getPad(2).getAxisZ().setLog(true);
myCanvas.getCanvas("genElectron").getPad(9).getAxisZ().setLog(true);


myCanvas.getCanvas("genKp1").divide(2,2);
myCanvas.getCanvas("genKp1").setGridX(false);
myCanvas.getCanvas("genKp1").setGridY(false);
myCanvas.getCanvas("genKp1").setAxisFontSize(18);
myCanvas.getCanvas("genKp1").setAxisTitleSize(24);
myCanvas.getCanvas("genKp1").draw(dg_kp1);
myCanvas.getCanvas("genKp1").getPad(0).getAxisZ().setLog(true);
myCanvas.getCanvas("genKp1").getPad(1).getAxisZ().setLog(true);


myCanvas.getCanvas("genKp2").divide(2,2);
myCanvas.getCanvas("genKp2").setGridX(false);
myCanvas.getCanvas("genKp2").setGridY(false);
myCanvas.getCanvas("genKp2").setAxisFontSize(18);
myCanvas.getCanvas("genKp2").setAxisTitleSize(24);
myCanvas.getCanvas("genKp2").draw(dg_kp2);
myCanvas.getCanvas("genKp2").getPad(0).getAxisZ().setLog(true);
myCanvas.getCanvas("genKp2").getPad(1).getAxisZ().setLog(true);

myCanvas.getCanvas("genKm").divide(2,2);
myCanvas.getCanvas("genKm").setGridX(false);
myCanvas.getCanvas("genKm").setGridY(false);
myCanvas.getCanvas("genKm").setAxisFontSize(18);
myCanvas.getCanvas("genKm").setAxisTitleSize(24);
myCanvas.getCanvas("genKm").draw(dg_km);
myCanvas.getCanvas("genKm").getPad(0).getAxisZ().setLog(true);
myCanvas.getCanvas("genKm").getPad(1).getAxisZ().setLog(true);

myCanvas.getCanvas("genProton").divide(2,2);
myCanvas.getCanvas("genProton").setGridX(false);
myCanvas.getCanvas("genProton").setGridY(false);
myCanvas.getCanvas("genProton").setAxisFontSize(18);
myCanvas.getCanvas("genProton").setAxisTitleSize(24);
myCanvas.getCanvas("genProton").draw(dg_pro);
myCanvas.getCanvas("genProton").getPad(0).getAxisZ().setLog(true);
myCanvas.getCanvas("genProton").getPad(1).getAxisZ().setLog(true);

myCanvas.getCanvas("genPim").divide(2,2);
myCanvas.getCanvas("genPim").setGridX(false);
myCanvas.getCanvas("genPim").setGridY(false);
myCanvas.getCanvas("genPim").setAxisFontSize(18);
myCanvas.getCanvas("genPim").setAxisTitleSize(24);
myCanvas.getCanvas("genPim").draw(dg_pim);
myCanvas.getCanvas("genPim").getPad(0).getAxisZ().setLog(true);
myCanvas.getCanvas("genPim").getPad(1).getAxisZ().setLog(true);

myCanvas.getCanvas("MM").divide(2,2);
myCanvas.getCanvas("MM").setGridX(false);
myCanvas.getCanvas("MM").setGridY(false);
myCanvas.getCanvas("MM").setAxisFontSize(18);
myCanvas.getCanvas("MM").setAxisTitleSize(24);
myCanvas.getCanvas("MM").draw(dg_mm);

JFrame frame = new JFrame("MC Study");
frame.setSize(1600, 1000);
frame.add(myCanvas);
frame.setLocationRelativeTo(null);
frame.setVisible(true);

	

	void fitxi(H1F hixi, F1D f1xi){

		// get histogram maximum in the rane 1.8 -1.86
		int i1 = hixi.getXaxis().getBin(1.82);
		int i2 = hixi.getXaxis().getBin(1.85);
		double hiMax = 0;
		int imax = i1;
		for(int i = i1; i <= i2; i2++){
			if (hiMax < hixi.getBinContent(i)){
				imax = i;
				hiMax = hixi.getBinContent(i);
			}
		}
		double mean = hixi.getDataX(imax); //hiw.getDataX(hiw.getMaximumBin());
        double amp  = hiMax;//hiw.getBinContent(hiw.getMaximumBin());
        double sigma = 0.05;
        f1xi.setParameter(0, amp);
        f1xi.setParameter(1, mean);
        f1xi.setParameter(2, sigma);
        double rmax = mean + 1.0 * Math.abs(sigma);
        double rmin = mean - 1.0 * Math.abs(sigma);
        f1xi.setRange(rmin, rmax);
        DataFitter.fit(f1xi, hixi, "Q"); //No options uses error for sigma 
        hixi.setFunction(null);
        mean = f1xi.getParameter(1);
        sigma = f1xi.getParameter(2);
        rmax = mean + 3.0 * Math.abs(sigma);
        rmin = mean - 3.0 * Math.abs(sigma);
        f1xi.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1xi, hixi, "Q"); //No options uses error for sigma 
        hixi.setFunction(null);

	}

