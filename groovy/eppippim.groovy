import java.io.*;
import java.util.*;
import org.jlab.io.hipo.*;
import java.text.SimpleDateFormat;
import org.jlab.jnp.physics.*;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Node;
//import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.pdg.PDGParticle;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.clas.physics.Vector3;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.ui.TBrowser;
import org.jlab.groot.data.*;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.LatexText;
import org.jlab.groot.ui.LatexTextTools;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
//import org.jlab.io.hipo.*;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import javax.swing.JFrame;
/*
import from hayward_coatjava_extensions 
import extended_kinematic_fitters.*; 
*/
/**
 * @author akhanal
 *
 */

 public class eppippim {

    public float  Eb, Mp;
    public LorentzVector VB, VT, Ve, VGS, VhadronSystm;
    public float STT;

    public float e_xB, e_Q2, e_W, e_virphoton;
    public List<Particle> fdels;
    public List<Particle> pips, pims, kps, kms, prots;
    public List<Particle> fpips, fpims, fkps, fkms, fprots; // = new ArrayList<Particle>();
    public List<Particle> cpips, cpims, ckps, ckms, cprots;

    public boolean found_eFD;

    public int e_ft_part_ind, pip_part_ind, pim_part_ind, kp_part_ind, km_part_ind, prot_part_ind;

    public H1F H_FT_W, H_FT_Q2, H_virphoton;
    public H2F H_FT_W_Q2;

    public H1F hi_deltap0, hi_deltatheta0, hi_deltaphi0;
    public H1F hi_deltap, hi_deltatheta, hi_deltaphi, hi_mm2, hi_mm2_corrected, hi_mm2_miss_cpip, hi_missing_e_corrected, hi_missing_e, hi_miss_e_cpip;

    public H2F hi_deltap_p, hi_deltap_theta, hi_deltap_phi, hi_deltap_phi1, hi_deltap_phi2, hi_deltap_phi3;
    public H2F hi_deltatheta_p, hi_deltatheta_theta, hi_deltatheta_phi, hi_deltatheta_phi1, hi_deltatheta_phi2, hi_deltatheta_phi3;
    public H2F hi_deltaphi_p, hi_deltaphi_theta, hi_deltaphi_phi, hi_deltaphi_phi1, hi_deltaphi_phi2, hi_deltaphi_phi3;


    public F1D fn_deltaphi_phi1, fn_deltaphi_phi2, fn_deltaphi_phi3, fn_deltaphi_theta, fn_deltaphi_p;
    public F1D fn_deltatheta_phi1, fn_deltatheta_phi2, fn_deltatheta_phi3, fn_deltatheta_theta, fn_deltatheta_p;
    public F1D fn_deltap_phi1, fn_deltap_phi2, fn_deltap_phi3, fn_deltap_theta, fn_deltap_p;
    public F1D fnp_cpip;

    public F1D fn2_deltaphi_phi1, fn2_deltaphi_phi2, fn2_deltaphi_phi3;
    public F1D fn2_deltatheta_theta, fn2_deltatheta_p;

    //Pierre's MC driven momentum correction for proton in FD (reg1:: < 27 deg, reg2:: >27deg, and, CD)
    public F1D fn_dpvsp_fdpro_reg1, fn_dpvsp_fdpro_reg2, fn_dpvsp_cdpro; 

    public EmbeddedCanvasTabbed myCanvas;

    public DataGroup dg_kinematics, dg_all,dg_dpdthdphi, dg_dp, dg_dth, dg_dphi;

    final int RED = 2;
    final int BLUE = 9;
    final int LIGHTGREEN = 3;
    final int LIGHTBROWN = 45;
    final int PINK = 46;

    public eppippim(float eb){
        Eb = eb;
        Mp = (float) PDGDatabase.getParticleById(2212).mass();
        VB = new LorentzVector(0, 0, Eb, Eb);
        VT = new LorentzVector(0, 0, 0, Mp);

        H_FT_W_Q2 = new H2F("H-FT-W-Q2", "H-FT-W-Q2", 100, 0, 1, 100, 0, 12);
        H_FT_W_Q2.setTitle("FT Q^2 vs W");
        H_FT_W_Q2.setTitleX("W ( GeV)");
        H_FT_W_Q2.setTitleY("Q^2 (GeV^2)");

        H_FT_W = new H1F("H-FT-W", "H-FT-W", 100, 0, 1);
        H_FT_W.setTitle("electron W");
        H_FT_W.setTitleX("W (GeV)");
        H_FT_W.setTitleY("count");
        H_FT_W.setFillColor(LIGHTGREEN);


        H_FT_Q2 = new H1F("H-FT-Q2", "H-FT-Q2", 100, 0, 12);
        H_FT_Q2.setFillColor(LIGHTGREEN);
        H_FT_Q2.setTitleX("Q^2 (GeV^2)");

        H_virphoton = new H1F("H-virphoton", "H-virphoton", 100, 0, 12);
        H_virphoton.setFillColor(LIGHTGREEN);
        H_virphoton.setTitleX("E_#gamma (GeV)");

        dg_kinematics = new DataGroup(2,2);
        dg_kinematics.addDataSet(H_FT_Q2, 0);
        dg_kinematics.addDataSet(H_FT_W, 1);
        dg_kinematics.addDataSet(H_FT_W_Q2, 2);
    }

    public eppippim(){
    //  NFTElec = 0;
    //  Eb = 10.575f;
    //  Eb = 10.604f;
        Eb = 10.1998f
    //  Eb = 7.54626f;
    //  Eb = 6.535f
        Mp = (float) PDGDatabase.getParticleById(2212).mass();
    //  Mp = 0.93827f;
        
        VB = new LorentzVector(0, 0, Eb, Eb);
        VT = new LorentzVector(0, 0, 0, Mp);


        fn_dpvsp_fdpro_reg1 = new F1D("fn_dpvsp_fdpro_reg1", "[a]+[b]*x+[c]*x*x", 0, 1);
        fn_dpvsp_fdpro_reg1.setParameters(0.153319, -0.298968, 0.1607);

        fn_dpvsp_fdpro_reg2 = new F1D("fn_dpvsp_fdpro_reg2", "[a]+[b]*x+[c]*x*x", 0, 1);
        fn_dpvsp_fdpro_reg2.setParameters(0.0398946,-0.0748125,0.0395764);

        fn_dpvsp_cdpro = new F1D("fn_dpvsp_cdpro", "[a]+[b]*x", 0, 1);
        fn_dpvsp_cdpro.setParameters(0.0292947,-0.0577956);


        H_FT_W_Q2 = new H2F("H-FT-W-Q2", "H-FT-W-Q2", 100, 0, 5, 100, 0, 12);
        H_FT_W_Q2.setTitle("FT Q^2 vs W");
        H_FT_W_Q2.setTitleX("W ( GeV)");
        H_FT_W_Q2.setTitleY("Q^2 (GeV^2)");
        
        
        H_FT_W = new H1F("H-FT-W", "H-FT-W", 100, 0, 5);
        H_FT_W.setTitle("electron W");
        H_FT_W.setTitleX("W (GeV)");
        H_FT_W.setTitleY("count");
        H_FT_W.setFillColor(LIGHTGREEN);


        H_FT_Q2 = new H1F("H-FT-Q2", "H-FT-Q2", 100, 0, 12);
        H_FT_Q2.setFillColor(LIGHTGREEN);
        H_FT_Q2.setTitleX("Q^2 (GeV^2)");

        H_virphoton = new H1F("H-virphoton", "H-virphoton", 100, 0, 12);
        H_virphoton.setFillColor(LIGHTGREEN);
        H_virphoton.setTitleX("E_#gamma (GeV)");

        dg_kinematics = new DataGroup(2,2);
        dg_kinematics.addDataSet(H_FT_Q2, 0);
        dg_kinematics.addDataSet(H_FT_W, 1);
        dg_kinematics.addDataSet(H_FT_W_Q2, 2);

        hi_mm2 = new H1F("hi-mm2", "hi-mm2", 100, -0.05, 0.05);
        hi_mm2.setTitleX("MM2(ep(f)#pi^-(f)#pi^+(c))");
        hi_mm2.setLineColor(RED);
        //hi_mm2.setFillColor(LIGHTGREEN);

        hi_missing_e = new H1F("hi-missing-e", "hi-missing-e", 100, -0.6, 0.6);
        hi_missing_e.setTitleX("missing E(ep(f)#pi^-(f)#pi^+(c))");
        hi_missing_e.setLineColor(RED);
        //hi_missing_e.setFillColor(LIGHTGREEN);

        hi_mm2_corrected = new H1F("hi-mm2-corrected", "hi-mm2-corrected", 100, -0.05, 0.05);
        hi_mm2_corrected.setTitleX("MM2(ep(f)#pi^-(f)#pi^+(c))");
        //hi_mm2_corrected.setFillColor(LIGHTGREEN);

        hi_missing_e_corrected = new H1F("hi-missing-e-corrected", "hi-missing-e-corrected", 100, -0.6, 0.6);
        hi_missing_e_corrected.setTitleX("missing E(ep(f)#pi^-(f)#pi^+(c))");
        //hi_missing_e_corrected.setFillColor(LIGHTGREEN);

        dg_all = new DataGroup(2,2);
        dg_all.addDataSet(hi_mm2, 0);
        dg_all.addDataSet(hi_mm2_corrected, 0);
        dg_all.addDataSet(hi_missing_e, 1);
        dg_all.addDataSet(hi_missing_e_corrected, 1);
        //dg_all.addDataSet(hi_mm2_corrected, 2);
        //dg_all.addDataSet(hi_missing_e_corrected, 2);
        
        

        hi_deltap0 = new H1F("hi-deltap0", "hi-deltap0", 100, -0.6, 0.6);
        hi_deltap0.setTitleX("#Deltap/p");
        //hi_deltap0.setFillColor(LIGHTGREEN);
        hi_deltap0.setLineColor(RED);

        hi_deltatheta0 = new H1F("hi-deltatheta0", "hi-deltatheta0", 100, -20, 20);
        hi_deltatheta0.setTitleX("#Delta#theta");
        //hi_deltatheta0.setFillColor(LIGHTGREEN);
        hi_deltatheta0.setLineColor(RED);

        hi_deltaphi0 = new H1F("hi-deltaphi0", "hi-deltaphi0", 100, -10, 10);
        hi_deltaphi0.setTitleX("#Delta#phi");
        //hi_deltaphi0.setFillColor(LIGHTGREEN);
        hi_deltaphi0.setLineColor(RED);


        hi_deltap = new H1F("hi-deltap", "hi-deltap", 100, -0.6, 0.6);
        hi_deltap.setTitleX("#Deltap/p");
        //hi_deltap.setFillColor(LIGHTGREEN);

        hi_deltatheta = new H1F("hi-deltatheta", "hi-deltatheta", 100, -20, 20); //100, -20, 20
        hi_deltatheta.setTitleX("#Delta#theta");
       //hi_deltatheta.setFillColor(LIGHTGREEN);

        hi_deltaphi = new H1F("hi-deltaphi", "hi-deltaphi", 100, -10, 10);
        hi_deltaphi.setTitleX("#Delta#phi");
        //hi_deltaphi.setFillColor(LIGHTGREEN);

        hi_mm2_miss_cpip = new H1F("hi-mm2-miss-cpip", "hi-mm2-miss-cpip", 100, -1, 1);
        hi_mm2_miss_cpip.setTitleX("MM2(efpf#pi^-)");
        hi_mm2_miss_cpip.setFillColor(LIGHTGREEN);

        hi_miss_e_cpip = new H1F("hi-miss-e-cpip", "hi-miss-e-cpip", 100, -1, 1);
        hi_miss_e_cpip.setTitleX("missing E(efpf#pi^-)")
        hi_miss_e_cpip.setFillColor(LIGHTGREEN);

        dg_dpdthdphi = new DataGroup(3,2);
        dg_dpdthdphi.addDataSet(hi_deltap, 0);
        dg_dpdthdphi.addDataSet(hi_deltap0, 0);
        dg_dpdthdphi.addDataSet(hi_deltatheta, 1);
        dg_dpdthdphi.addDataSet(hi_deltatheta0, 1);
        dg_dpdthdphi.addDataSet(hi_deltaphi, 2);
        dg_dpdthdphi.addDataSet(hi_deltaphi0, 2);
        dg_dpdthdphi.addDataSet(hi_mm2_miss_cpip, 3);
        dg_dpdthdphi.addDataSet(hi_miss_e_cpip, 4);

        hi_deltap_p = new H2F("hi-deltap-p", "hi-deltap-p", 4, 0.3, 1.1, 30, -0.15, 0.15);//30, 0.3, 1.2,
        hi_deltap_p.setTitleX("p");
        hi_deltap_p.setTitleY("#Deltap/p");
        hi_deltap_theta = new H2F("hi-deltap-theta", "hi-deltap-theta", 6, 35, 60, 30, -0.15, 0.15);//40, 35, 70,//5, 35, 60
        hi_deltap_theta.setTitleX("#theta");
        hi_deltap_theta.setTitleY("#Deltap/p");
        hi_deltap_phi = new H2F("hi-deltap-phi", "hi-deltap-phi", 12, -180, 180, 30, -0.15, 0.15);//60, -180, 180//6, -180, 180
        hi_deltap_phi.setTitleX("#phi");
        hi_deltap_phi.setTitleY("#Deltap/p");
        hi_deltap_phi1 = new H2F("hi-deltap-phi1", "hi-deltap-phi1", 12, -180, 180, 30, -0.15, 0.15);
        hi_deltap_phi1.setTitleX("#phi1 + 90^o");
        hi_deltap_phi1.setTitleY("#Deltap/p");
        hi_deltap_phi2 = new H2F("hi-deltap-phi2", "hi-deltap-phi2", 12, -180, 180, 30, -0.15, 0.15);
        hi_deltap_phi2.setTitleX("#phi2");
        hi_deltap_phi2.setTitleY("#Deltap/p");
        hi_deltap_phi3 = new H2F("hi-deltap-phi3", "hi-deltap-phi3", 12, -180, 180, 30, -0.15, 0.15);
        hi_deltap_phi3.setTitleX("#phi3");
        hi_deltap_phi3.setTitleY("#Deltap/p");
        /*
        fn_deltap_phi1 = new F1D("fn-deltap-phi1", "[a]+[b]*x+[c]*x*x", -150, 0);
        fn_deltap_phi1.setParameter(0, 0.03791121);
        fn_deltap_phi1.setParameter(1, 0.001047737);
        fn_deltap_phi1.setParameter(2, 0.000006543885);
        fn_deltap_phi1.setOptStat("1111");
        fn_deltap_phi1.setLineColor(RED);
        */
        fn_deltap_phi1 = new F1D("fn-deltap-phi1", "[a]+[b]*x", -150, 0);
        fn_deltap_phi1.setParameters(-0.01437840, 0.0001040396);
        /*
        fn_deltap_phi2 = new F1D("fn-deltap-phi2", "[a]+[b]*x", -110, 20);
        fn_deltap_phi2.setParameter(0, 0.03685392);
        fn_deltap_phi2.setParameter(1, -0.0006716383);
        fn_deltap_phi2.setOptStat("1111");
        fn_deltap_phi2.setLineColor(RED);
        */
        fn_deltap_phi2 = new F1D("fn-deltap-phi2", "[a]+[b]*x", -110, 20);
        fn_deltap_phi2.setParameters(0.01956305, -0.0006708716);
        /*
        fn_deltap_phi3 = new F1D("fn-deltap-phi3", "[a]+[b]*x+[c]*x*x", 10, 150);
        fn_deltap_phi3.setParameter(0, 0.08084716);
        fn_deltap_phi3.setParameter(1, -0.001602984);
        fn_deltap_phi3.setParameter(2, 0.000007449761);
        fn_deltap_phi3.setOptStat("11111");
        fn_deltap_phi3.setLineColor(RED);
        */
        fn_deltap_phi3 = new F1D("fn-deltap-phi3", "[a]+[b]*x", 10, 150);
        fn_deltap_phi3.setParameters(0.02556854, -0.0005140661);
        /*
        fn_deltap_theta = new F1D("fn-deltap-theta", "[a]+[b]*x+[c]*x*x", 35, 70);
        fn_deltap_theta.setParameter(0, 0.1257472);//0.09477161
        fn_deltap_theta.setParameter(1, -0.005097940);//-0.003668270
        fn_deltap_theta.setParameter(2, 0.00005068624);//0.00003604987
        fn_deltap_theta.setOptStat("1111");
        fn_deltap_theta.setLineColor(RED);
        */
        fn_deltap_theta = new F1D("fn-deltap-theta", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 35, 70);
        fn_deltap_theta.setParameters(0.3052909,  -0.01665354, 0.0002970831, -0.000001718886);//0.09477161

        /*
        fn_deltap_p = new F1D("fn-deltap-p", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 0.3, 1.2);
        fn_deltap_p.setParameter(0, -0.0008726564);
        fn_deltap_p.setParameter(1, -0.01093897);
        fn_deltap_p.setParameter(2, 0.04739136);
        fn_deltap_p.setParameter(3, -0.02852069);
        fn_deltap_p.setOptStat("11111");
        fn_deltap_p.setLineColor(RED);
        */
        fn_deltap_p = new F1D("fn-deltap-p", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 0.3, 1.2);
        fn_deltap_p.setParameters(0.06089941, -0.2664393, 0.3686122, -0.1617709);

        dg_dp = new DataGroup(6, 2);
        dg_dp.addDataSet(hi_deltap_p, 0);
        dg_dp.addDataSet(hi_deltap_theta, 1);
        dg_dp.addDataSet(hi_deltap_phi, 2);
        dg_dp.addDataSet(hi_deltap_phi1, 3);
        dg_dp.addDataSet(hi_deltap_phi2, 4);
        dg_dp.addDataSet(hi_deltap_phi3, 5);
        //dg_dp.addDataSet(fn_deltap_p, 6);
        //dg_dp.addDataSet(fn_deltap_theta, 7);
        //dg_dp.addDataSet(fn_deltap_p, 8);
        //dg_dp.addDataSet(fn_deltap_phi1, 9);
        //dg_dp.addDataSet(fn_deltap_phi2, 10);
        //dg_dp.addDataSet(fn_deltap_phi3, 11);



       //hi_deltatheta_p = new H2F("hi-deltatheta-p", "hi-deltatheta-p", 30, 0.3, 1.2, 30, -70, 70);
        //hi_deltatheta_p = new H2F("hi-deltatheta-p", "hi-deltatheta-p", 30, 0.3, 1.2, 30, -10, 10);
        hi_deltatheta_p = new H2F("hi-deltatheta-p", "hi-deltatheta-p", 4, 0.3, 1.1, 30, -10, 10); //30, 0.3, 1.2, 30, -10, 10
        hi_deltatheta_p.setTitleX("p");
        hi_deltatheta_p.setTitleY("#Delta#theta");
        hi_deltatheta_theta = new H2F("hi-deltatheta-theta", "hi-deltatheta-theta", 6, 35, 60, 30, -10, 10);
        hi_deltatheta_theta.setTitleX("#theta");
        hi_deltatheta_theta.setTitleY("#Delta#theta");
        hi_deltatheta_phi = new H2F("hi-deltatheta-phi", "hi-deltatheta-phi", 12, -180, 180, 30, -10, 10);
        hi_deltatheta_phi.setTitleX("#phi1");
        hi_deltatheta_phi.setTitleY("#Delta#theta");
        hi_deltatheta_phi1 = new H2F("hi-deltatheta-phi1", "hi-deltatheta-phi1", 12, -180, 180, 30, -10, 10);
        hi_deltatheta_phi1.setTitleX("#phi1 + 90^o");
        hi_deltatheta_phi1.setTitleY("#Delta#theta");
        hi_deltatheta_phi2 = new H2F("hi-deltatheta-phi2", "hi-deltatheta-phi2", 12, -180, 180, 30, -10, 10);
        hi_deltatheta_phi2.setTitleX("#phi2");
        hi_deltatheta_phi2.setTitleY("#Delta#theta");
        hi_deltatheta_phi3 = new H2F("hi-deltatheta-phi3", "hi-deltatheta-phi3", 12, -180, 180, 30, -10, 10);
        hi_deltatheta_phi3.setTitleX("#phi3");
        hi_deltatheta_phi3.setTitleY("#Delta#theta");

        // theta correction functions for three regions of CVT from deltatheta Vs phi plot
        /*
        fn_deltatheta_phi1 = new F1D("fn-deltatheta-phi1", "[a]+[b]*x", -150, 0);
        fn_deltatheta_phi1.setParameter(0, 0.4797968);
        fn_deltatheta_phi1.setParameter(1, 0.0047468);
        fn_deltatheta_phi1.setOptStat("1111");
        fn_deltatheta_phi1.setLineColor(RED);
        */
        fn_deltatheta_phi1 = new F1D("fn-deltatheta-phi1", "[a]+[b]*x", -150, 0);
        fn_deltatheta_phi1.setParameters(-0.2754198, 0.008523673);

        /*
        fn_deltatheta_phi2 = new F1D("fn-deltatheta-phi2", "[a]+[b]*x", -110, 20);
        fn_deltatheta_phi2.setParameter(0, -0.1036288);
        fn_deltatheta_phi2.setParameter(1, -0.0089518);
        fn_deltatheta_phi2.setOptStat("1111");
        fn_deltatheta_phi2.setLineColor(RED);
        */
        fn_deltatheta_phi2 = new F1D("fn-deltatheta-phi2", "[a]+[b]*x", -110, 20);
        fn_deltatheta_phi2.setParameters(-1.112222, -0.01297936);

        /*
        fn_deltatheta_phi3 = new F1D("fn-deltatheta-phi3", "[a]+[b]*x", 10, 150);
        fn_deltatheta_phi3.setParameter(0, 0.071482);
        fn_deltatheta_phi3.setParameter(1, 0.0002006);
        fn_deltatheta_phi3.setOptStat("1111");
        fn_deltatheta_phi3.setLineColor(RED);
        */
        fn_deltatheta_phi3 = new F1D("fn-deltatheta-phi3", "[a]+[b]*x", 10, 150);
        fn_deltatheta_phi3.setParameters(-0.7841433, -0.002032877);

         /*
        fn_deltatheta_theta = new F1D("fn-deltatheta-theta", "[a]+[b]*x", 35, 70);//37, 65
        fn_deltatheta_theta.setParameter(0, 1.402174);
        fn_deltatheta_theta.setParameter(1, -0.07413168);
        //*/
        /*
        fn_deltatheta_theta = new F1D("fn-deltatheta-theta", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 35, 70);//37, 65
        fn_deltatheta_theta.setParameter(0, 31.61607);//31.61607 //31.1886
        fn_deltatheta_theta.setParameter(1, -1.843130);//-1.843130 //-1.810522
        fn_deltatheta_theta.setParameter(2, 0.03357832);//0.03357832 //0.0329588
        fn_deltatheta_theta.setParameter(3, -0.0002073033);//-0.0002073033 //-0.000203
        fn_deltatheta_theta.setOptStat("11111");
        fn_deltatheta_theta.setLineColor(RED);
       // */

        fn_deltatheta_theta = new F1D("fn-deltatheta-theta", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 35, 70);//37, 65
        fn_deltatheta_theta.setParameters(35.99945, -2.163462, 0.04280008, -0.0002816622);//31.61607 //31.1886

        fn2_deltatheta_theta = new F1D("fn2-deltatheta-theta", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 35, 70);//37, 65
        fn2_deltatheta_theta.setParameters(14.05754, -0.8545750, 0.01715881, -0.0001160209);//31.61607 //31.1886

        /*
        fn_deltatheta_p = new F1D("fn-deltatheta-p", "[a]+[b]*x", 0.3, 1.2);
        fn_deltatheta_p.setParameter(0, 1.106319);
        fn_deltatheta_p.setParameter(1, -1.326232);
        fn_deltatheta_p.setOptStat("1111");
        fn_deltatheta_p.setLineColor(RED);
        */
        fn_deltatheta_p = new F1D("fn-deltatheta-p", "[a]+[b]*x", 0.3, 1.2);
        fn_deltatheta_p.setParameters(0.8548455, -0.7965661);
        fn2_deltatheta_p = new F1D("fn2-deltatheta-p", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 0.3, 1.2);
        fn2_deltatheta_p.setParameters(2.229385, -7.408962, 9.294125, -3.932364);
        

        dg_dth = new DataGroup(6, 2);
        dg_dth.addDataSet(hi_deltatheta_p, 0);
        dg_dth.addDataSet(hi_deltatheta_theta, 1);
        dg_dth.addDataSet(hi_deltatheta_phi, 2);
        dg_dth.addDataSet(hi_deltatheta_phi1, 3);
        dg_dth.addDataSet(hi_deltatheta_phi2, 4);
        dg_dth.addDataSet(hi_deltatheta_phi3, 5);
        //dg_dth.addDataSet(fn_deltatheta_p, 6);
        //dg_dth.addDataSet(fn_deltatheta_theta, 7);
        //dg_dth.addDataSet(hi_deltatheta_phi, 8);
        //dg_dth.addDataSet(fn_deltatheta_phi1, 9);
        //dg_dth.addDataSet(fn_deltatheta_phi2, 10);
        //dg_dth.addDataSet(fn_deltatheta_phi3, 11);

        hi_deltaphi_p = new H2F("hi-deltaphi-p", "hi-deltaphi-p", 4, 0.3, 1.1, 30, -4, 4);//
        hi_deltaphi_p.setTitleX("p");
        hi_deltaphi_p.setTitleY("#Delta#phi");
        hi_deltaphi_theta = new H2F("hi-deltaphi-theta", "hi-deltaphi-theta", 6, 35, 60, 30, -4, 4);
        hi_deltaphi_theta.setTitleX("#theta");
        hi_deltaphi_theta.setTitleY("#Delta#phi");
        hi_deltaphi_phi = new H2F("hi-deltaphi-phi", "hi-deltaphi-phi", 12, -180, 180, 30, -4, 4);
        hi_deltaphi_phi.setTitleX("#phi");
        hi_deltaphi_phi.setTitleY("#Delta#phi");
        hi_deltaphi_phi1 = new H2F("hi-deltaphi-phi1", "hi-deltaphi-phi1", 12, -180, 180, 30, -4, 4);
        hi_deltaphi_phi1.setTitleX("#phi1 + 90^o");
        hi_deltaphi_phi1.setTitleY("#Delta#phi");
        hi_deltaphi_phi2 = new H2F("hi-deltaphi-phi2", "hi-deltaphi-phi2", 12, -180, 180, 30, -4, 4);
        hi_deltaphi_phi2.setTitleX("#phi2");
        hi_deltaphi_phi2.setTitleY("#Delta#phi");
        hi_deltaphi_phi3 = new H2F("hi-deltaphi-phi3", "hi-deltaphi-phi3", 12, -180, 180, 30, -4, 4);
        hi_deltaphi_phi3.setTitleX("#phi3");
        hi_deltaphi_phi3.setTitleY("#Delta#phi");

        // phi correction functions for three regions of CVT from deltaPhi Vs phi plot
        /*
        fn_deltaphi_phi1 = new F1D("fn-deltaphi-phi1", "[a]+[b]*x", -150, 0);
        fn_deltaphi_phi1.setParameter(0, 0.3560393);//0.3560393
        fn_deltaphi_phi1.setParameter(1, 0.001261624);//0.001261624
        fn_deltaphi_phi1.setOptStat("1111");
        fn_deltaphi_phi1.setLineColor(RED);
        */
        fn_deltaphi_phi1 = new F1D("fn-deltaphi-phi1", "[a]+[b]*x", -150, 0);
        fn_deltaphi_phi1.setParameters(0.642, 0.005);//0.3560393
        fn2_deltaphi_phi1 = new F1D("fn2-deltaphi-phi1", "[a]+[b]*x", -150, 0);
        fn2_deltaphi_phi1.setParameters(-0.04600986, -0.0004717446);

        /*
        fn_deltaphi_phi2 = new F1D("fn-deltaphi-phi2", "[a]+[b]*x", -110, 20);
        fn_deltaphi_phi2.setParameter(0, -0.8892721);
        fn_deltaphi_phi2.setParameter(1, -0.01833468);
        fn_deltaphi_phi2.setOptStat("1111");
        fn_deltaphi_phi2.setLineColor(RED);
        */
        fn_deltaphi_phi2 = new F1D("fn-deltaphi-phi2", "[a]+[b]*x", -110, 20);
        fn_deltaphi_phi2.setParameters(-0.986, -0.019);
        fn2_deltaphi_phi2 = new F1D("fn2-deltaphi-phi2", "[a]+[b]*x", -110, 20);
        fn2_deltaphi_phi2.setParameters(-0.1833572, -0.002527996);

        /*
        fn_deltaphi_phi3 = new F1D("fn-deltaphi-phi3", "[a]+[b]*x", 10, 150);
        fn_deltaphi_phi3.setParameter(0, -1.007084);
        fn_deltaphi_phi3.setParameter(1, 0.005664863);
        fn_deltaphi_phi3.setOptStat("1111");
        fn_deltaphi_phi3.setLineColor(RED);
        */
        fn_deltaphi_phi3 = new F1D("fn-deltaphi-phi3", "[a]+[b]*x", 10, 150);
        fn_deltaphi_phi3.setParameters(-0.987, 0.006);

        fn2_deltaphi_phi3 = new F1D("fn2-deltaphi-phi3", "[a]+[b]*x", 10, 150);
        fn2_deltaphi_phi3.setParameters(-0.09040472, -0.0005510582);


       /*
        fnp_cpip = new F1D("fnp_cpip", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 0.3, 1.2);
        fnp_cpip.setParameter(0, 0.03105106);
        fnp_cpip.setParameter(1, -0.05107082);
        fnp_cpip.setParameter(2, 0.1000333);
        fnp_cpip.setParameter(3, -0.05812612);

       //*/

       fn_deltaphi_theta = new F1D("fn_deltaphi_theta", "[a]+[b]*x+[c]*x*x+[d]*x*x*x", 35, 60);
       fn_deltaphi_theta.setParameters(9.355429, -0.5630870, 0.01123901, -0.0000742);

       fn_deltaphi_p = new F1D("fn_deltaphi_p", "[a]+[b]*x+[c]*x*x", 0.3, 1.1);
       fn_deltaphi_p.setParameters(1.190272, -2.473123, 1.261091);

        

        dg_dphi = new DataGroup(6, 2);
        dg_dphi.addDataSet(hi_deltaphi_p, 0);
        dg_dphi.addDataSet(hi_deltaphi_theta, 1);
        dg_dphi.addDataSet(hi_deltaphi_phi, 2);
        dg_dphi.addDataSet(hi_deltaphi_phi1, 3);
        dg_dphi.addDataSet(hi_deltaphi_phi2, 4);
        dg_dphi.addDataSet(hi_deltaphi_phi3, 5);
        //dg_dphi.addDataSet(hi_deltaphi_p, 0);
        //dg_dphi.addDataSet(hi_deltaphi_theta, 1);
        //dg_dphi.addDataSet(hi_deltaphi_phi, 2);
        //dg_dphi.addDataSet(fn_deltaphi_phi1, 9);
        //dg_dphi.addDataSet(fn_deltaphi_phi2, 10);
        //dg_dphi.addDataSet(fn_deltaphi_phi3, 11);



    }

    public void fillRecBank(DataBank recBank) {
        STT = recBank.getFloat("startTime", 0);
    }

    public int makeFDElectron(DataBank bank) {

        int NFDElec = 0;
        fdels = new ArrayList<Particle>();
        Particle Vfdel = null;
        for(int k = 0; k < bank.rows(); k++){
            int pid = bank.getInt("pid", k);
            byte q = bank.getByte("charge", k);
            float px = bank.getFloat("px", k);
            float py = bank.getFloat("py", k);
            float pz = bank.getFloat("pz", k);
            float vx = bank.getFloat("vx", k);
            float vy = bank.getFloat("vy", k);
            float vz = bank.getFloat("vz", k);
            int partstatus = bank.getShort("status", k);
            if (partstatus<0) partstatus = -partstatus;
            boolean inDC = (partstatus >= 2000 && partstatus < 4000);
            float e_mom = (float) Math.sqrt(px * px + py * py + pz * pz);
            float e_the = (float) Math.toDegrees(Math.acos(pz / e_mom));
            if (pid == 11 && inDC && bank.getShort("status", k) < 0.0 && q == -1){
                NFDElec++;
                Vfdel = new Particle(pid, px, py, pz, vx, vy, vz);
                fdels.add(Vfdel);
            }
            if (NFDElec == 1) {
                found_eFD = true;
                double eft = Math.sqrt(e_mom * e_mom + PDGDatabase.getParticleById(11).mass()*PDGDatabase.getParticleById(11).mass());
                float e_phi = (float) Math.toDegrees(Math.atan2(py, px));
                Ve = new LorentzVector(px, py, pz, eft);
                VGS = new LorentzVector(0, 0, 0, 0);
                VGS.add(VB);
                VGS.sub(Ve);
                e_Q2 = (float) -VGS.mass2();
                VhadronSystm = new LorentzVector(0, 0, 0, 0);
                VhadronSystm.add(VB);
                VhadronSystm.add(VT);
                VhadronSystm.sub(Ve);
                e_xB = e_Q2 / (2f * Mp * (Eb - e_mom));
                e_W = (float) VhadronSystm.mass();
                e_virphoton = (float) (Eb - Ve.e());
                //e_W = (float) Math.sqrt(Mp * Mp + e_Q2 * (1f / e_xB - 1f));
                return k;

            }

        }
        return -1;

    }
    public void makeOthers(DataBank recbank, DataBank recSCBank) {

        int nkp = 0;
        int npip = 0;
        int nfdpip = 0;
        int nfdwithcutpip = 0;
        int ncdpip = 0;
        int ncdwithcutpip = 0;
        
        // particles going central + forward
        pips = new ArrayList<Particle>();
        pims = new ArrayList<Particle>();
        kps = new ArrayList<Particle>();
        kms = new ArrayList<Particle>();
        prots = new ArrayList<Particle>();
        
        // particles going forward
        fpips = new ArrayList<Particle>();
        fpims = new ArrayList<Particle>();
        fkps = new ArrayList<Particle>();
        fkms = new ArrayList<Particle>();
        fprots = new ArrayList<Particle>();
        
        // particles going central
        cpips = new ArrayList<Particle>();
        cpims = new ArrayList<Particle>();
        ckps = new ArrayList<Particle>();
        ckms = new ArrayList<Particle>();
        cprots = new ArrayList<Particle>();
        
        for (int k = 0; k < recbank.rows(); k++) {
            byte q = recbank.getByte("charge", k);
            int pid = recbank.getInt("pid", k);
            //int ftbpid = recFTbank.getInt("pid", k);
            float px = recbank.getFloat("px", k);
            float py = recbank.getFloat("py", k);
            float pz = recbank.getFloat("pz", k);
            float vx = recbank.getFloat("vx", k);
            float vy = recbank.getFloat("vy", k);
            float vz = recbank.getFloat("vz", k);
            float chi2pid = recbank.getFloat("chi2pid", k);
            float beta = recbank.getFloat("beta", k);
            //float ftbbe = recFTbank.getFloat("beta", k);
            int status = recbank.getShort("status", k);
            if (status<0) status = -status;
            boolean inDC = (status >= 2000 && status < 4000);
            boolean inCD = (status >= 4000);
            float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
            float the = (float) Math.toDegrees(Math.acos(pz / mom));
            float phi = (float) Math.toDegrees(Math.atan2(py, px));
            float mass = mom * mom * ((1 -  beta * beta)/ ( beta * beta ));
            //float FTBmass = mom * mom * (1 / ( ftbbe * ftbbe ) - 1);
            //boolean fdkpMass2Cut = (0.2 < FTBmass && FTBmass < 0.40); // cut for K- candidates in FD 
            //boolean cdkpMass2Cut = (0.15 < FTBmass && FTBmass < 0.35); //  cut for K- candidates in CD 
            //boolean kpMass2Cut = (0.180625 < FTBmass && FTBmass < 0.36);// 0.425 < kpM < 0.6
            boolean fdChi2pidCut = (Math.abs(chi2pid) < 5.0);
            boolean cdChi2pidCut = (Math.abs(chi2pid) < 5.0);
            if (pid == 211 ) {
                npip++;
                pip_part_ind = k;
                
                for (int r = 0; r < recSCBank.rows(); r++) {
                    if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {   
                        if (recSCBank.getShort("pindex", r) == pip_part_ind && recSCBank.getByte("layer", r) == 2) {    
                            nfdpip++;
                            //pip_FTOF_pad1b = recSCBank.getShort("component", r);
                            float pip_FTOF1b_t = recSCBank.getFloat("time", r);
                            //pip_FTOF1b_path = recSCBank.getFloat("path", r);
                            float pip_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass());
                            //pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_ftb_beta * 29.98f) - STT - pip_vz/ (pip_ftb_beta * 29.98f);
                            //pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_beta * 29.98f) - STT - pip_vz/ (pip_beta * 29.98f);
                            //pip_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
                            float pip_FTOF1b_TOF = (float) pip_FTOF1b_t - STT - vz/ (pip_beta * 29.98f);
                            //difference in the measured and computed vertex times 
                            float pip_FTOF1b_vt = (float) pip_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass())));
                            if (Math.abs(pip_FTOF1b_vt) < 0.5) { //0.4
                                nfdwithcutpip++;
                                Particle Vpip = new Particle(pid, px, py, pz, vx, vy, vz);
                                fpips.add(Vpip);
                                pips.add(Vpip);
                                //H_pip_vt_p.fill(pip_FTOF1b_vt, pip_mom);
                            }
                            /*
                            H_pip_vt_p.fill(pip_FTOF1b_vt, pip_mom);
                            H_pip_vz.fill(pip_vz);
                            H_pip_FTOF1b_t.fill(pip_FTOF1b_TOF);
                            H_pip_FTOF1b_path.fill(pip_FTOF1b_path);
                            */
                            
                        } // pip from FTOF panal 1b
                    }
                    
                    if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut) {
                        
                        if (recSCBank.getShort("pindex", r) == pip_part_ind ) {
                            ncdpip++;
                            //pip_CTOF_pad = recSCBank.getShort("component", r);
                            float pip_CTOF_t = recSCBank.getFloat("time", r);
                            //pip_CTOF_path = recSCBank.getFloat("path", r);
                            float pipc_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass());
                            //pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_ftb_beta * 29.98f) - STT - pip_vz/ (pip_ftb_beta * 29.98f);
                            //pip_CTOF_vt = pip_CTOF_t - pip_CTOF_path / (pipc_beta * 29.98f) - STT - pip_vz/ (pipc_beta * 29.98f);
                            //pip_CTOF_vt = STT - recFTbank.getFloat("vt", k);
                            float pip_CTOF_TOF = (float) pip_CTOF_t - STT - vz/ (pipc_beta * 29.98f);
                            //difference in the measured and computed vertex times 
                            float pip_CTOF_vt = (float) pip_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass())));
                            if (Math.abs(pip_CTOF_vt) < 0.4) {
                                ncdwithcutpip++;
                                Particle Vpipc = new Particle(pid, px, py, pz, vx, vy, vz);  
                                pips.add(Vpipc);
                                cpips.add(Vpipc);
                                //H_pip_vt_p.fill(pip_FTOF1b_vt, pip_mom);
                            }

                            /*
                            H_pipc_vt_p.fill(pip_CTOF_vt, pip_mom);
                            H_pipc_vz.fill(pip_vz);
                            H_pip_CTOF_t.fill(pip_CTOF_t - STT - pip_vz/ (pipc_beta * 29.98f));
                            H_pip_CTOF_path.fill(pip_CTOF_path);
                            */
                            
                            
                        }
                        
                        
                    } //CTOF
                    
                    
                }           

            }
            if (pid == -211 ) {
                pim_part_ind = k;
                
                for (int r = 0; r < recSCBank.rows(); r++) {
                    if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {   
                        if (recSCBank.getShort("pindex", r) == pim_part_ind && recSCBank.getByte("layer", r) == 2) {
                            //pim_FTOF_pad1b = recSCBank.getShort("component", r);
                            float pim_FTOF1b_t = recSCBank.getFloat("time", r);
                            //pim_FTOF1b_path = recSCBank.getFloat("path", r);
                            float pim_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass());
                            //pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_ftb_beta * 29.98f) - STT - pim_vz/ (pim_ftb_beta * 29.98f);
                            //pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_beta * 29.98f) - STT - pim_vz/ (pim_beta * 29.98f);
                            //pim_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
                            float pim_FTOF1b_TOF = (float) pim_FTOF1b_t - STT - vz/ (pim_beta * 29.98f);
                            //difference in the measured and computed vertex times 
                            float pim_FTOF1b_vt = (float) pim_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass())));
                            if (Math.abs(pim_FTOF1b_vt) < 0.5) {
                                //Vpim = new LorentzVector(pim_px*pimCor, pim_py*pimCor, pim_pz*pimCor, Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass()));
                                Particle Vpim = new Particle(pid, px, py, pz, vx, vy, vz);
                                pims.add(Vpim);
                                fpims.add(Vpim);
                                //H_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);
                            }
                            /*
                            H_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);
                            H_pim_vz.fill(pim_vz);
                            H_pim_FTOF1b_t.fill(pim_FTOF1b_t - STT - pim_vz/ (pim_beta * 29.98f));
                            H_pim_FTOF1b_path.fill(pim_FTOF1b_path);
                            */
                        } // pim from FTOF panal 1b
                    }
                    
                    if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut ) {
                        
                        if (recSCBank.getShort("pindex", r) == pim_part_ind ) {
                            //pim_CTOF_pad = recSCBank.getShort("component", r);
                            float pim_CTOF_t = recSCBank.getFloat("time", r);
                            //pim_CTOF_path = recSCBank.getFloat("path", r);
                            float pimc_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass());
                            //pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_ftb_beta * 29.98f) - STT - pim_vz/ (pim_ftb_beta * 29.98f);
                            //pim_CTOF_vt = pim_CTOF_t - pim_CTOF_path / (pimc_beta * 29.98f) - STT - pim_vz/ (pimc_beta * 29.98f);
                            //pim_CTOF_vt = STT - recFTbank.getFloat("vt", k);
                            float pim_CTOF_TOF = (float) pim_CTOF_t - STT - vz/ (pimc_beta * 29.98f);
                            //difference in the measured and computed vertex times 
                            float pim_CTOF_vt = (float) pim_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass())));
                            if (Math.abs(pim_CTOF_vt) < 0.4) {
                                //Vpimc = new LorentzVector(pim_px*pimCor, pim_py*pimCor, pim_pz*pimCor, Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass()));
                                Particle Vpimc = new Particle(pid, px, py, pz, vx, vy, vz);      
                                pims.add(Vpimc);
                                cpims.add(Vpimc);
                                //H_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);
                            }
                            /*
                            H_pimc_vt_p.fill(pim_CTOF_vt, pim_mom);
                            H_pimc_vz.fill(pim_vz);
                            H_pim_CTOF_t.fill(pim_CTOF_t - STT - pim_vz/ (pimc_beta * 29.98f));
                            H_pim_CTOF_path.fill(pim_CTOF_path);
                            */
                            
                            
                        }
                        
                        
                    } //CTOF
                    
                    
                }
            }
            if (pid == 321 ) {
                kp_part_ind = k;
                
                for (int r = 0; r < recSCBank.rows(); r++) {
                    if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut ) {  //&& fdkpMass2Cut
                        if (recSCBank.getShort("pindex", r) == kp_part_ind && recSCBank.getByte("layer", r) == 2) {
                            //kp_FTOF_pad1b = recSCBank.getShort("component", r);
                            float kp_FTOF1b_t = recSCBank.getFloat("time", r);
                            //kp_FTOF1b_path = recSCBank.getFloat("path", r);
                            float kp_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass());
                            //kp_FTOF1b_vt = kp_FTOF1b_t - kp_FTOF1b_path / (k_ftb_beta * 29.98f) - STT - kp_vz/ (kp_ftb_beta * 29.98f);
                            //kp_FTOF1b_vt = kp_FTOF1b_t - kp_FTOF1b_path / (kp_beta * 29.98f) - STT - kp_vz/ (kp_beta * 29.98f);
                            //kp_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
                            float kp_FTOF1b_TOF = (float) kp_FTOF1b_t - STT - vz/ (kp_beta * 29.98f);
                            //difference in the measured and computed vertex times 
                            float kp_FTOF1b_vt = (float) kp_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass())));
                            boolean fdkpvzCut = vz > -10 && vz < 2;
                            if (Math.abs(kp_FTOF1b_vt) < 0.5 ) {// && kp_mom < 2.8 //fdkpvzCut && 
                                Particle Vkp = new Particle(pid, px, py, pz, vx, vy, vz);
                                fkps.add(Vkp);
                                kps.add(Vkp);
                                //H_kp_vt_p.fill(kp_FTOF1b_vt, kp_mom);
                            }
                            /*
                            H_kp_vt_p.fill(kp_FTOF1b_vt, kp_mom);
                            H_kp_vz.fill(kp_vz);
                            H_kp_FTOF1b_t.fill(kp_FTOF1b_t - STT - kp_vz/ (kp_beta * 29.98f));
                            H_kp_FTOF1b_path.fill(kp_FTOF1b_path);
                            */
                        } // kp from FTOF panal 1b
                    }
                    
                    if (recSCBank.getByte("detector", r) == 4 && inCD  && cdChi2pidCut ) { //&& cdkpMass2Cut
                        
                        if (recSCBank.getShort("pindex", r) == kp_part_ind ) {
                            //kp_CTOF_pad = recSCBank.getShort("component", r);
                            float kp_CTOF_t = recSCBank.getFloat("time", r);
                            //kp_CTOF_path = recSCBank.getFloat("path", r);
                            float kpc_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass());
                            //kp_CTOF_vt = kp_CTOF_t - kp_CTOF_path / (kpc_beta * 29.98f) - STT - kp_vz/ (kpc_beta * 29.98f);
                            //kp_CTOF_vt = STT - recFTbank.getFloat("vt", k);
                            boolean cdkpvzCut = vz > -8 && vz < 1;
                            float kp_CTOF_TOF = (float) kp_CTOF_t - STT - vz/ (kpc_beta * 29.98f);
                            //difference in the measured and computed vertex times 
                            float kp_CTOF_vt = (float) kp_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass())));
                            if (Math.abs(kp_CTOF_vt) < 0.4) { //cdkpvzCut &&
                                Particle Vkpc = new Particle(pid, px, py, pz, vx, vy, vz);
                                kps.add(Vkpc);
                                ckps.add(Vkpc);
                                //H_kp_vt_p.fill(kp_FTOF1b_vt, kp_mom);
                            }
                            /*
                            H_kpc_vt_p.fill(kp_CTOF_vt, kp_mom);
                            H_kpc_vz.fill(kp_vz);
                            H_kp_CTOF_t.fill(kp_CTOF_t - STT - kp_vz/ (kpc_beta * 29.98f));
                            H_kp_CTOF_path.fill(kp_CTOF_path);
                            */
                            
                            
                        }
                        
                        
                    } //CTOF
                    
                    
                }
             //   nkp++;

            }
            
            if (pid == -321) {
                km_part_ind = k;
                
                for (int r = 0; r < recSCBank.rows(); r++) {
                    if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {   
                        if (recSCBank.getShort("pindex", r) == km_part_ind && recSCBank.getByte("layer", r) == 2) {
                            //km_FTOF_pad1b = recSCBank.getShort("component", r);
                            float km_FTOF1b_t = recSCBank.getFloat("time", r);
                            //km_FTOF1b_path = recSCBank.getFloat("path", r);
                            float km_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass());
                            //km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_ftb_beta * 29.98f) - STT - km_vz/ (km_ftb_beta * 29.98f);
                            //km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_beta * 29.98f) - STT - km_vz/ (km_beta * 29.98f);
                            //km_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
                            boolean fdkmvzCut = vz > -10 && vz < 2;
                            float km_FTOF1b_TOF = (float) km_FTOF1b_t - STT - vz/ (km_beta * 29.98f);
                            boolean km_FTOF1b_TOFCut = 22 < km_FTOF1b_TOF && km_FTOF1b_TOF < 28;
                            //difference in the measured and computed vertex times 
                            float km_FTOF1b_vt = (float) km_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass())));

                            if (Math.abs(km_FTOF1b_vt) < 0.5) {//&& fdkmvzCut && km_FTOF1b_TOFCut
                                Particle Vkm = new Particle(pid, px, py, pz, vx, vy, vz);
                                fkms.add(Vkm);
                                kms.add(Vkm);
                                //H_km_vt_p.fill(km_FTOF1b_vt, km_mom);
                            }
                            /*
                            H_km_vt_p.fill(km_FTOF1b_vt, km_mom);
                            H_km_vz.fill(km_vz);
                            H_km_FTOF1b_t.fill(km_FTOF1b_TOF);
                            H_km_FTOF1b_path.fill(km_FTOF1b_path);
                            */
                        } // km from FTOF panal 1b
                    }
                    
                    if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut ) {
                        
                        if (recSCBank.getShort("pindex", r) == km_part_ind ) {
                            //km_CTOF_pad = recSCBank.getShort("component", r);
                            float km_CTOF_t = recSCBank.getFloat("time", r);
                            //km_CTOF_path = recSCBank.getFloat("path", r);
                            float kmc_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass());
                            //km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_ftb_beta * 29.98f) - STT - km_vz/ (km_ftb_beta * 29.98f);
                            //km_CTOF_vt = km_CTOF_t - km_CTOF_path / (kmc_beta * 29.98f) - STT - km_vz/ (kmc_beta * 29.98f);
                            //km_CTOF_vt = STT - recFTbank.getFloat("vt", k);
                            boolean cdkmvzCut = vz > -8 && vz < 1;
                            float km_CTOF_TOF = (float) km_CTOF_t - STT - vz/ (kmc_beta * 29.98f);
                            //difference in the measured and computed vertex times 
                            float km_CTOF_vt = (float) km_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass())));
                            if (Math.abs(km_CTOF_vt) < 0.4 ) { //&& cdkmvzCut
                                Particle Vkmc = new Particle(pid, px, py, pz, vx, vy, vz);     
                                kms.add(Vkmc);
                                ckms.add(Vkmc);
                                //H_km_vt_p.fill(km_FTOF1b_vt, km_mom);
                            }
                            /*
                            H_kmc_vt_p.fill(km_CTOF_vt, km_mom);
                            H_kmc_vz.fill(km_vz);
                            H_km_CTOF_t.fill(km_CTOF_t - STT - km_vz/ (kmc_beta * 29.98f));
                            H_km_CTOF_path.fill(km_CTOF_path);
                            */
                            
                            
                        }
                        
                        
                    } //CTOF                
                }
            }
            
              
             
            
            if (pid == 2212 ) {
                prot_part_ind = k;
                
                for (int r = 0; r < recSCBank.rows(); r++) {

                    if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {   
                        if (recSCBank.getShort("pindex", r) == prot_part_ind && recSCBank.getByte("layer", r) == 2) {
                            //prot_FTOF_pad1b = recSCBank.getShort("component", r);
                            float prot_FTOF1b_t = recSCBank.getFloat("time", r);
                            //prot_FTOF1b_path = recSCBank.getFloat("path", r);
                            float prot_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass());
                            //prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_ftb_beta * 29.98f) - STT - prot_vz/ (prot_ftb_beta * 29.98f);
                            //prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_beta * 29.98f) - STT - prot_vz/ (prot_beta * 29.98f);
                            //prot_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
                            float prot_FTOF1b_TOF = (float) prot_FTOF1b_t - STT - vz/ (prot_beta * 29.98f);
                            boolean prot_FTOF1b_TOFCut = 22 < prot_FTOF1b_TOF && prot_FTOF1b_TOF < 32;
                            //difference in the measured and computed vertex times 
                            float prot_FTOF1b_vt = (float) prot_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass())));
                            if (Math.abs(prot_FTOF1b_vt) < 0.5 ) {// && prot_FTOF1b_TOFCut
                                //Vprot = new LorentzVector(prot_px * proCor, prot_py* proCor, prot_pz* proCor, Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass()));
                                Particle Vprot = new Particle(pid, px, py, pz, vx, vy, vz);
                                //Pierre's MC driven correction 
                                if (Math.toDegrees(Vprot.theta()) > 27.0) Vprot.setP(Vprot.p() + fn_dpvsp_fdpro_reg1.evaluate(Vprot.p()));
                                if (Math.toDegrees(Vprot.theta()) < 27.0) Vprot.setP(Vprot.p() + fn_dpvsp_fdpro_reg2.evaluate(Vprot.p()));

                                fprots.add(Vprot);
                                prots.add(Vprot);
                                //H_prot_vt_p.fill(prot_FTOF1b_vt, prot_mom);
                            }
                            /*
                            H_prot_vt_p.fill(prot_FTOF1b_vt, prot_mom);
                            H_prot_vz.fill(prot_vz);
                            H_prot_FTOF1b_t.fill(prot_FTOF1b_t - STT - prot_vz/ (prot_beta * 29.98f));
                            H_prot_FTOF1b_path.fill(prot_FTOF1b_path);
                            */
                        } // prot from FTOF panal 1b
                    }
                    
                    if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut ) {
                        
                        if (recSCBank.getShort("pindex", r) == prot_part_ind ) {
                            //prot_CTOF_pad = recSCBank.getShort("component", r);
                            float prot_CTOF_t = recSCBank.getFloat("time", r);
                            //prot_CTOF_path = recSCBank.getFloat("path", r);
                            float protc_beta = mom / (float) Math.sqrt(mom * mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass());
                            //prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_ftb_beta * 29.98f) - STT - prot_vz/ (prot_ftb_beta * 29.98f);
                            //prot_CTOF_vt = prot_CTOF_t - prot_CTOF_path / (protc_beta * 29.98f) - STT - prot_vz/ (protc_beta * 29.98f);
                            //prot_CTOF_vt = STT - recFTbank.getFloat("vt", k);
                            float prot_CTOF_TOF = (float) prot_CTOF_t - STT - vz/ (protc_beta * 29.98f);
                            //difference in the measured and computed vertex times 
                            float prot_CTOF_vt = (float) prot_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass())));
                            if (Math.abs(prot_CTOF_vt) < 0.4) {
                                //Vprotc = new LorentzVector(prot_px * proCor, prot_py* proCor, prot_pz* proCor, Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass()));
                                Particle Vprotc = new Particle(pid, px, py, pz, vx, vy, vz); 
                                //Pierre's MC driven correction
                                Vprotc.setP(Vprotc.p() + fn_dpvsp_cdpro.evaluate(Vprotc.p()));    
                                cprots.add(Vprotc);
                                prots.add(Vprotc);
                                //H_prot_vt_p.fill(prot_FTOF1b_vt, prot_mom);
                            }
                            /*
                            H_protc_vt_p.fill(prot_CTOF_vt, prot_mom);
                            H_protc_vz.fill(prot_vz);
                            H_prot_CTOF_t.fill(prot_CTOF_t - STT - prot_vz/ (protc_beta * 29.98f));
                            H_prot_CTOF_path.fill(prot_CTOF_path);
                            */
                            
                        }
                        
                        
                    } //CTOF
                        
                }
            
            }

        }// For loop

    }// end Makeother

    public void processEvent(DataEvent event) {
        resetCounters();

        //if (event.hasBank("RECFT::Event"))
        //  fillRecBank(event.getBank("RECFT::Event"));
        if (event.hasBank("REC::Event"))
            fillRecBank(event.getBank("REC::Event"));
        //if (event.hasBank("REC::Particle") && event.hasBank("RECFT::Particle") && event.hasBank("REC::ForwardTagger")) e_ft_part_ind = makeFTElectron(event.getBank("REC::Particle"), event.getBank("RECFT::Particle"), event.getBank("REC::ForwardTagger"));
        if (event.hasBank("REC::Particle")) e_ft_part_ind = makeFDElectron(event.getBank("REC::Particle"));
        //if (event.hasBank("MC::Particle") == true) fillMCPartBank(event.getBank("MC::Particle"));
        if(e_ft_part_ind > -1 && found_eFD) {       
            if (event.hasBank("REC::Particle") && event.hasBank("REC::Scintillator")) makeOthers(event.getBank("REC::Particle"), event.getBank("REC::Scintillator"));
            
            //if (event.hasBank("REC::Scintillator")) fillFTOF(event.getBank("REC::Particle"), event.getBank("REC::Scintillator"));
            //System.out.println("found pip :: " + pips.size() + " pips tracks in FD & " + cpips.size() + " pip tracks in CD.");
            FillHists();
            
        } // e_ft_part_ind > -1
        

    }

    public void resetCounters() {

        found_eFD = false;
        e_ft_part_ind = -1;
        pip_part_ind = -1;
        pim_part_ind = -1;
        kp_part_ind = -1;
        km_part_ind = -1;
        prot_part_ind = -1;

       // /*
        fdels = null;
        kps = null;
        pips = null;
        kms = null;
        pims = null;
        prots = null;
        fkps = null;
        fpips = null;
        fkms = null;
        fpims = null;
        fprots = null;
        ckps = null;
        cpips = null;
        ckms = null;
        cpims = null;
        cprots = null;
     //   */

    }
    public void analyze() {

    }

    public boolean select_efpcpipfpim(){
        boolean res = false;
        if(found_eFD && cpips.size() == 1 && fprots.size() == 1 && fpims.size() == 1){
            if(Ve.e()>0.1*Eb && fprots.get(0).e()>0.94358 && fpims.get(0).e() > 0.3 && cpips.get(0).e() > 0.2){

                LorentzVector VRHO0 = new LorentzVector(0,0,0,0);
                VRHO0.add(cpips.get(0).vector());
                VRHO0.add(fpims.get(0).vector());

                LorentzVector Q = new LorentzVector(0,0,0,0);
                Q.add(VB);
                Q.sub(Ve);
                LorentzVector W = new LorentzVector(0,0,0,0);
                W.add(Q);
                W.add(VT);

                if( -Q.mass2()>0.8 && W.mass()>1.8 ){
                    LorentzVector VmissP = new LorentzVector(0,0,0,0);
                    VmissP.add(W);
                    VmissP.sub(VRHO0);

                    LorentzVector VmissPIM = new LorentzVector(0,0,0,0);
                    VmissPIM.add(W);
                    VmissPIM.sub(fprots.get(0).vector());
                    VmissPIM.sub(cpips.get(0).vector());

                    LorentzVector VmissPIP = new LorentzVector(0,0,0,0);
                    VmissPIP.add(W);
                    VmissPIP.sub(fprots.get(0).vector());
                    VmissPIP.sub(fpims.get(0).vector());

                    LorentzVector VmissAll = new LorentzVector(0,0,0,0);
                    VmissAll.add(VmissPIP);
                    VmissAll.sub(cpips.get(0).vector());




                    res = true 
                        && VmissAll.e() > -0.5 && VmissAll.e() < 0.5 
                        && VmissAll.mass2() > -0.1 &&  VmissAll.mass2() < 0.1 
                        && VmissAll.px()*VmissAll.px() + VmissAll.py()*VmissAll.py() < 0.5 
                        //&& Vangle( fprots.get(0).vector().vect() , VmissP.vect() ) < 12
                        && Vangle( cpips.get(0).vector().vect() , VmissPIP.vect() ) < 12
                        //&& Vangle( fpims.get(0).vector().vect() , VmissPIM.vect() ) < 12
                        && VmissP.mass2() > 0.5 && VmissP.mass2() < 1.5 
                        && VmissPIP.mass2() > -0.3 && VmissPIP.mass2() < 0.5 
                        && VmissPIM.mass2() > -0.3 && VmissPIM.mass2() < 0.5
                        ;
                } 
                //res = true;
            }
        }
        return res;
    }

    public double Vangle(Vector3 v1, Vector3 v2){ 
        double res=0;
        double l1 = v1.mag();
        double l2 = v2.mag();
        if( l1*l2 > 0)res = Math.toDegrees( Math.acos( v1.dot(v2)/(l1*l2) ) );
        return res;
    }

    
    public double cdpart_thetaCorrn(Particle cdpart){
        double newTheta = -1000;
        double newTheta1 = -1000;
        double newTheta2 = -1000;
        double oldTheta = Math.toDegrees(cdpart.theta());
        boolean correctPhi = true;
        boolean correctTheta = true;
        boolean correctP = true;

        //phi dependent theta corrections

        if(Math.toDegrees(cdpart.phi())>135 || Math.toDegrees(cdpart.phi()) < -105){
            double phi1 = Math.toDegrees(cdpart.phi());

            if(Math.toDegrees(cdpart.phi()) > 135){
                phi1 = Math.toDegrees(cdpart.phi()) - 270;

            }
            else{
                phi1 = Math.toDegrees(cdpart.phi()) + 90;
            }
            //newPhi = oldPhi * (1 - fn_deltaphi_phi1.evaluate(phi1));
            newTheta = oldTheta + fn_deltatheta_phi1.evaluate(phi1);
        }

        if(Math.toDegrees(cdpart.phi()) > -105 && Math.toDegrees(cdpart.phi()) < 15){
            double phi1 = Math.toDegrees(cdpart.phi());

            newTheta = oldTheta + fn_deltatheta_phi2.evaluate(phi1);
        }
        if(Math.toDegrees(cdpart.phi())>15 && Math.toDegrees(cdpart.phi()) < 135){
            double phi1 = Math.toDegrees(cdpart.phi());

            newTheta = oldTheta + fn_deltatheta_phi3.evaluate(phi1);
        }
        
                   
        if(correctTheta){//&& correctTheta
            cdpart.setTheta(Math.toRadians(newTheta));
        }

        // theta dependent theta correction
        if(Math.toDegrees(cdpart.theta())>0 ){//&& Math.toDegrees(cpips.get(0).theta()) < 70
            //newTheta1 = newTheta + fn_deltatheta_theta.evaluate(oldTheta);
            newTheta1 = newTheta + fn_deltatheta_theta.evaluate(newTheta);
        }
                    
     //   /*    
        if(correctTheta){
            //cpips.get(0).vector().vect().setMagThetaPhi(cpips.get(0).p(), Math.toRadians(newTheta), cpips.get(0).phi());
            cdpart.setTheta(Math.toRadians(newTheta1));
        }

        if(cdpart.p()>0 ){//&& Math.toDegrees(cpips.get(0).theta()) < 70
            newTheta2 = newTheta1 + fn_deltatheta_p.evaluate(cdpart.p());
        }
        /*
        if(correctTheta){
            //cpips.get(0).vector().vect().setMagThetaPhi(cpips.get(0).p(), Math.toRadians(newTheta), cpips.get(0).phi());
            cdpart.setTheta(Math.toRadians(newTheta2));
        }
        */
        return newTheta2;
    }

    public double cdpart_thetaCorrn_ita2(Particle cdpart){
        double newTheta = -1000;
        double newTheta1 = -1000;
        double newTheta2 = -1000;
        double oldTheta = Math.toDegrees(cdpart.theta());
        boolean correctPhi = true;
        boolean correctTheta = true;
        boolean correctP = true;

        // theta dependent theta correction
        if(Math.toDegrees(cdpart.theta())>0 ){//&& Math.toDegrees(cpips.get(0).theta()) < 70
            //newTheta1 = newTheta + fn_deltatheta_theta.evaluate(oldTheta);
            newTheta1 = oldTheta + fn2_deltatheta_theta.evaluate(oldTheta);
        }
                    
     //   /*    
        if(correctTheta){
            //cpips.get(0).vector().vect().setMagThetaPhi(cpips.get(0).p(), Math.toRadians(newTheta), cpips.get(0).phi());
            cdpart.setTheta(Math.toRadians(newTheta1));
        }

        
        if(cdpart.p()>0 ){//&& Math.toDegrees(cpips.get(0).theta()) < 70
            newTheta2 = newTheta1 + fn2_deltatheta_p.evaluate(cdpart.p());
        }

        return newTheta2;

    }


    public double cdpart_phiCorrn(Particle cdpart){

        boolean correctPhi = true;

        double newPhi = -1000;
        double newPhi1 = -1000;
        double newPhi2 = -1000;
        double oldPhi = Math.toDegrees(cdpart.phi());
        double oldTheta = Math.toDegrees(cdpart.theta());
        double oldP = cdpart.p();
        // phi dependent phi correction
        if(Math.toDegrees(cdpart.phi())>135 || Math.toDegrees(cdpart.phi()) < -105){

            double phi1 = Math.toDegrees(cdpart.phi());

                if(Math.toDegrees(cdpart.phi()) > 135){
                    phi1 = Math.toDegrees(cdpart.phi()) - 270;

                }
                else{
                    phi1 = Math.toDegrees(cdpart.phi()) + 90;
                }
                       //newPhi = oldPhi * (1 - fn_deltaphi_phi1.evaluate(phi1));
            newPhi = oldPhi + fn_deltaphi_phi1.evaluate(phi1);
        }

        if(Math.toDegrees(cdpart.phi()) > -105 && Math.toDegrees(cdpart.phi()) < 15){

            newPhi = oldPhi + fn_deltaphi_phi2.evaluate(oldPhi);
        }
        if(Math.toDegrees(cdpart.phi())>15 && Math.toDegrees(cdpart.phi()) < 135){
            newPhi = oldPhi + fn_deltaphi_phi3.evaluate(oldPhi);
        }
        // update particle phi
        if(correctPhi ){
            cdpart.vector().vect().setMagThetaPhi(cdpart.p(), cdpart.theta(), Math.toRadians(newPhi));
        }
        // theta dependent cpip phi correction
        if(correctPhi){
            newPhi1 = newPhi + fn_deltaphi_theta.evaluate(Math.toDegrees(cdpart.theta()));
        }
        // update particle phi
        if(correctPhi ){
            cdpart.vector().vect().setMagThetaPhi(cdpart.p(), cdpart.theta(), Math.toRadians(newPhi1));
        }

        if(correctPhi){
            newPhi2 = newPhi1 + fn_deltaphi_p.evaluate(cdpart.p());
        }

    /*
        if(correctPhi ){//&& correctTheta
            cpips.get(0).vector().vect().setMagThetaPhi(cpips.get(0).p(), cpips.get(0).theta(), Math.toRadians(newPhi));
        }
        */
        return newPhi2;


    }

    public double cdpart_phiCorrn_it2(Particle cdpart){
        boolean correctPhi = true;

        double newPhi = -1000;
        double newPhi1 = -1000;
        double newPhi2 = -1000;
        double oldPhi = Math.toDegrees(cdpart.phi());
        double oldTheta = Math.toDegrees(cdpart.theta());
        double oldP = cdpart.p();
        // phi dependent phi correction
        if(Math.toDegrees(cdpart.phi())>135 || Math.toDegrees(cdpart.phi()) < -105){

            double phi1 = Math.toDegrees(cdpart.phi());

                if(Math.toDegrees(cdpart.phi()) > 135){
                    phi1 = Math.toDegrees(cdpart.phi()) - 270;

                }
                else{
                    phi1 = Math.toDegrees(cdpart.phi()) + 90;
                }
                       //newPhi = oldPhi * (1 - fn_deltaphi_phi1.evaluate(phi1));
            newPhi = oldPhi + fn2_deltaphi_phi1.evaluate(phi1);
        }

        if(Math.toDegrees(cdpart.phi()) > -105 && Math.toDegrees(cdpart.phi()) < 15){

            newPhi = oldPhi + fn2_deltaphi_phi2.evaluate(oldPhi);
        }
        if(Math.toDegrees(cdpart.phi())>15 && Math.toDegrees(cdpart.phi()) < 135){
            newPhi = oldPhi + fn2_deltaphi_phi3.evaluate(oldPhi);
        }
        // update particle phi
        if(correctPhi ){
            cdpart.vector().vect().setMagThetaPhi(cdpart.p(), cdpart.theta(), Math.toRadians(newPhi));
        }
        return newPhi;


    }

    public double cdpart_pCorrn(Particle cdpart){
        boolean correctPhi = true;
        boolean correctTheta = true;
        boolean correctP = true;
        // phi dependent P correction
        double newP = -1000;
        double newP1 = -1000;
        double newP2 = -1000;
        double oldP = cdpart.p();

        if(Math.toDegrees(cdpart.phi())>135 || Math.toDegrees(cdpart.phi()) < -105){

            double phi1 = Math.toDegrees(cdpart.phi());

            if(Math.toDegrees(cdpart.phi()) > 135){
                phi1 = Math.toDegrees(cdpart.phi()) - 270;

            }
            else{
                phi1 = Math.toDegrees(cdpart.phi()) + 90;
            }
            newP = oldP*(1 + fn_deltap_phi1.evaluate(phi1));
        }

        if(Math.toDegrees(cdpart.phi()) > -105 && Math.toDegrees(cdpart.phi()) < 15){
            double phi1 = Math.toDegrees(cdpart.phi());

            newP = oldP*(1 + fn_deltap_phi2.evaluate(phi1));
        }
        if(Math.toDegrees(cdpart.phi())>15 && Math.toDegrees(cdpart.phi()) < 135){
            double phi1 = Math.toDegrees(cdpart.phi());

            newP = oldP*(1 + fn_deltap_phi3.evaluate(phi1));
        }
                   
        if(correctP){
            cdpart.setP(newP);
        }

        // theta dependent cpip p correction
        if(Math.toDegrees(cdpart.theta())>0){
            newP1 = newP*(1 + fn_deltap_theta.evaluate(Math.toDegrees(cdpart.theta())));
        }
        if(correctP){
            cdpart.setP(newP1);
        }
        /*
        // p dependent cpip p correction
        
        if(cdpart.p()>0.0){
            newP2 = newP1*(1 + fn_deltap_p.evaluate(newP1));
        }
        if(correctP){
            cdpart.setP(newP2);
        }
        
        /*

        if(correctP){
            cdpart.setP(newP1);

        } 
        */
        return  newP1;             
    }

    public void FillHists() {

        if(found_eFD){      
            H_FT_W_Q2.fill(e_W, e_Q2);
            H_FT_W.fill(e_W);
            H_FT_Q2.fill(e_Q2); 
            H_virphoton.fill(e_virphoton);  
        }

        if(select_efpcpipfpim()){

            LorentzVector VRHO0 = new LorentzVector(0,0,0,0);
            VRHO0.add(cpips.get(0).vector());
            VRHO0.add(fpims.get(0).vector());

            LorentzVector Q = new LorentzVector(0,0,0,0);
            Q.add(VB);
            Q.sub(Ve);
            LorentzVector W = new LorentzVector(0,0,0,0);
            W.add(Q);
            W.add(VT);

            if( -Q.mass2()>0.8 && W.mass()>1.5 ){
                LorentzVector VmissCPIP = new LorentzVector(0,0,0,0);
                VmissCPIP.add(W);
                VmissCPIP.sub(fprots.get(0).vector());
                VmissCPIP.sub(fpims.get(0).vector());

                LorentzVector VmissAll = new LorentzVector(0,0,0,0);
                VmissAll.add(VmissCPIP);
                VmissAll.sub(cpips.get(0).vector());

                hi_mm2.fill(VmissAll.mass2());
                hi_missing_e.fill(VmissAll.e());

                if(VmissAll.mass2() > -0.1 && VmissAll.mass2() < 0.1 && VmissAll.e() > -0.5 && VmissAll.e() < 0.5){
                    hi_mm2_miss_cpip.fill(VmissCPIP.mass2());
                    hi_miss_e_cpip.fill(VmissCPIP.e()-cpips.get(0).e());
                    hi_deltap0.fill((VmissCPIP.p()-cpips.get(0).p())/cpips.get(0).p());
                    hi_deltatheta0.fill(Math.toDegrees(VmissCPIP.theta())-Math.toDegrees(cpips.get(0).theta()));
                    hi_deltaphi0.fill(Math.toDegrees(VmissCPIP.phi())-Math.toDegrees(cpips.get(0).phi()));

                   // System.out.println("cpipBeforeCorrn------- p = "+ cpips.get(0).p()+" theta = "+ cpips.get(0).theta() + " phi = "+ cpips.get(0).phi());

                    // appply corrections to phi first and theta correcotion 
                //   /*
                    boolean correctPhi = false;
                    boolean correctTheta = false;
                    boolean correctP = false;
                    //*/
                   /*
                    boolean correctPhi = true;
                    boolean correctTheta = true;
                    boolean correctP = true;
                   // */
                 //  /*
                    boolean correctPhi_ita2 = false;
                    boolean correctTheta_ita2 = false;
                    //*/

                   /*
                    boolean correctPhi_ita2 = true;
                    boolean correctTheta_ita2 = true;
                    //*/
                    if(correctPhi ){//&& correctTheta

                        cpips.get(0).vector().vect().setMagThetaPhi(cpips.get(0).p(), cpips.get(0).theta(), Math.toRadians(cdpart_phiCorrn(cpips.get(0))));
                        //cpips.get(0).vector().vect().setMagThetaPhi(cpips.get(0).p(), cpips.get(0).theta(), Math.toRadians(newPhi));
                    }

                    if(correctTheta){
                        //cpips.get(0).vector().vect().setMagThetaPhi(cpips.get(0).p(), Math.toRadians(newTheta), cpips.get(0).phi());
                        cpips.get(0).setTheta(Math.toRadians(cdpart_thetaCorrn(cpips.get(0))));
                    }

                    if(correctP){
                        cpips.get(0).setP(cdpart_pCorrn(cpips.get(0)));

                    }

                    if(correctPhi_ita2){
                        cpips.get(0).vector().vect().setMagThetaPhi(cpips.get(0).p(), cpips.get(0).theta(), Math.toRadians(cdpart_phiCorrn_it2(cpips.get(0))));
                    }

                    if(correctTheta_ita2){
                        cpips.get(0).setTheta(Math.toRadians(cdpart_thetaCorrn_ita2(cpips.get(0))));
                    }

                   // System.out.println("cpipBeforeCorrn------- p = "+ cpips.get(0).p()+" theta = "+ cpips.get(0).theta() + " phi = "+ cpips.get(0).phi());

                    float deltap = (VmissCPIP.p()-cpips.get(0).p())/cpips.get(0).p(); 
                    //float pres = deltap/cpips.get(0).p(); 
                    float deltatheta = Math.toDegrees(VmissCPIP.theta())-Math.toDegrees(cpips.get(0).theta());
                    float deltaphi = Math.toDegrees(VmissCPIP.phi())-Math.toDegrees(cpips.get(0).phi());

                    hi_deltap.fill(deltap);
                    hi_deltatheta.fill(deltatheta);
                    hi_deltaphi.fill(deltaphi);

                    hi_deltap_p.fill(cpips.get(0).p(), deltap);
                    hi_deltap_theta.fill(Math.toDegrees(cpips.get(0).theta()), deltap);
                    hi_deltap_phi.fill(Math.toDegrees(cpips.get(0).phi()), deltap);
                    

                    hi_deltatheta_p.fill(cpips.get(0).p(), deltatheta);
                    hi_deltatheta_theta.fill(Math.toDegrees(cpips.get(0).theta()), deltatheta);
                    hi_deltatheta_phi.fill(Math.toDegrees(cpips.get(0).phi()), deltatheta);

                    hi_deltaphi_p.fill(cpips.get(0).p(), deltaphi);
                    hi_deltaphi_theta.fill(Math.toDegrees(cpips.get(0).theta()), deltaphi);
                    hi_deltaphi_phi.fill(Math.toDegrees(cpips.get(0).phi()), deltaphi);

                    if(Math.toDegrees(cpips.get(0).phi())>135 || Math.toDegrees(cpips.get(0).phi()) < -105){

                        double phi1 = Math.toDegrees(cpips.get(0).phi());

                     //   /*
                        if(Math.toDegrees(cpips.get(0).phi()) > 135){
                            phi1 = Math.toDegrees(cpips.get(0).phi()) - 270;
                        }
                        else{
                            phi1 = Math.toDegrees(cpips.get(0).phi()) + 90;
                        }
                        //*/

                        hi_deltap_phi1.fill(phi1, deltap);
                        hi_deltatheta_phi1.fill(phi1, deltatheta);
                        hi_deltaphi_phi1.fill(phi1, deltaphi);
                    }

                    if(Math.toDegrees(cpips.get(0).phi()) > -105 && Math.toDegrees(cpips.get(0).phi()) < 15){
                        hi_deltap_phi2.fill(Math.toDegrees(cpips.get(0).phi()), deltap);
                        hi_deltatheta_phi2.fill(Math.toDegrees(cpips.get(0).phi()), deltatheta);
                        hi_deltaphi_phi2.fill(Math.toDegrees(cpips.get(0).phi()), deltaphi);

                    }
                    if(Math.toDegrees(cpips.get(0).phi())>15 && Math.toDegrees(cpips.get(0).phi()) < 135){
                        hi_deltap_phi3.fill(Math.toDegrees(cpips.get(0).phi()), deltap);
                        hi_deltatheta_phi3.fill(Math.toDegrees(cpips.get(0).phi()), deltatheta);
                        hi_deltaphi_phi3.fill(Math.toDegrees(cpips.get(0).phi()), deltaphi);


                    }
                    //*/

                    // check MM2 and Missing energy after correction
                    LorentzVector VmissAllCorrect = new LorentzVector(0,0,0,0);
                    VmissAllCorrect.add(VmissCPIP);
                    VmissAllCorrect.sub(cpips.get(0).vector());

                    hi_mm2_corrected.fill(VmissAllCorrect.mass2());
                    hi_missing_e_corrected.fill(VmissAllCorrect.e());


                }
                

            }

        }

    }


    public void sliceFit(H2F h2, F1D f1, GraphErrors meanX){
        meanX.reset();
        ArrayList<H1F> hslice = h2.getSlicesX();
        for(int i=0; i<hslice.size(); i++) {
            double  x = h2.getXAxis().getBinCenter(i);
            double ex = 0;
            double  y = hslice.get(i).getRMS();
            double ey = 0;
            double mean  = hslice.get(i).getDataX(hslice.get(i).getMaximumBin());
            double amp   = hslice.get(i).getBinContent(hslice.get(i).getMaximumBin());
            double sigma = hslice.get(i).getRMS();
            f1.setParameter(0, amp);
            f1.setParameter(1, mean);
            f1.setParameter(2, sigma);
            DataFitter.fit(f1, hslice.get(i), "Q"); //No options uses error for sigma 
            if(amp>0) meanX.addPoint(x, f1.getParameter(1), ex, f1.parameter(1).error());

        }
    }

    public void dataPointsBetweenSlices(H2F h2, F1D f1, ArrayList<GraphErrors> dataPoints){


        
        ArrayList<H1F> hslice = h2.getSlicesX();
        ArrayList<double> x = new ArrayList<double>(); //(hslice.size())???
        ArrayList<double> y = new ArrayList<double>();
        ArrayList<double> ex = new ArrayList<double>();
        ArrayList<double> ey = new ArrayList<double>();

        for(int i=0; i<hslice.size(); i++) {
           
            double mean  = hslice.get(i).getDataX(hslice.get(i).getMaximumBin());
            double amp   = hslice.get(i).getBinContent(hslice.get(i).getMaximumBin());
            double sigma = hslice.get(i).getRMS();
            f1.setParameter(0, amp);
            f1.setParameter(1, mean);
            f1.setParameter(2, sigma);
            DataFitter.fit(f1, hslice.get(i), "Q"); //No options uses error for sigma 
            if(amp>0){
                x.add(h2.getXAxis().getBinCenter(i));
                double errx = 0;
                ex.add(errx);//ex.add(h2.getXAxis().getBinWidth(i));
                y.add(f1.getParameter(1));
                ey.add(f1.parameter(1).error());

            }
            
        }

        //GraphErrors gr;
        //gr.reset();
        //dataPoints.reset();
        //dataPoints = new ArrayList<GraphErrors>;
        for(int i=0; i<(x.size()-1); i++){
            GraphErrors gr = new GraphErrors();
            gr.addPoint(x[i], y[i], ex[i], ey[i]);
            gr.addPoint(x[i+1], y[i+1], ex[i+1], ey[i+1])

            dataPoints.add(gr);
            //gr.reset();

        }



    }


    public void plotGraph(){
        // for tracks reconstructed twice in FD and CD
        H2F h1 = dg_dp.getH2F("hi-deltap-p");
        //H2F h1 = hi_fkpckp_deltap_p;
        F1D f1 = new F1D("f1","[amp]*gaus(x,[mean],[sigma])", -0.150, 0.150);
        GraphErrors gr1 = new GraphErrors();
        sliceFit(h1, f1, gr1);
        //DataFitter.fit(fn_deltap_p, gr1, "Q");
        ArrayList<GraphErrors> graphs = new ArrayList<GraphErrors>();
        dataPointsBetweenSlices(h1, f1, graphs);
        myCanvas.getCanvas("DeltaP").cd(6);
        //myCanvas.getCanvas("DeltaP").draw(gr1);
        myCanvas.getCanvas("DeltaP").draw(graphs[0]);
        //myCanvas.getCanvas("DeltaP").cd(6);
        for(int i=1; i< graphs.size(); i++){
            myCanvas.getCanvas("DeltaP").draw(graphs[i], "same");
        }
        //myCanvas.getCanvas("DeltaP").draw(graphs[graphs.size()-1], "same");
        //DataFitter.fit(fn_deltap_p, gr1, "Q");
        //myCanvas.getCanvas("DeltaP").cd(6);
        //myCanvas.getCanvas("DeltaP").draw(fn_deltap_p);

        H2F h2 = dg_dp.getH2F("hi-deltap-theta");
        //H2F h2 = hi_fkpckp_deltap_theta;
        F1D f2 = new F1D("f2","[amp]*gaus(x,[mean],[sigma])", -0.150, 0.150);
        GraphErrors gr2 = new GraphErrors();
        sliceFit(h2, f2, gr2);
        //DataFitter.fit(fn_deltap_theta, gr2, "Q");
        myCanvas.getCanvas("DeltaP").cd(7);
        myCanvas.getCanvas("DeltaP").draw(gr2);

        H2F h3 = dg_dp.getH2F("hi-deltap-phi");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f3 = new F1D("f3","[amp]*gaus(x,[mean],[sigma])", -0.150, 0.150);
        GraphErrors gr3 = new GraphErrors();
        sliceFit(h3, f3, gr3);
        myCanvas.getCanvas("DeltaP").cd(8);
        myCanvas.getCanvas("DeltaP").draw(gr3);

        H2F h31 = dg_dp.getH2F("hi-deltap-phi1");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f31 = new F1D("f31","[amp]*gaus(x,[mean],[sigma])", -0.150, 0.150);
        GraphErrors gr31 = new GraphErrors();
        sliceFit(h31, f31, gr31);
        //DataFitter.fit(fn_deltap_phi1, gr31, "Q");
        //fn_deltap_phi1.setParameter(0, fn_deltap_phi1.getParameter(0));
        myCanvas.getCanvas("DeltaP").cd(9);
        myCanvas.getCanvas("DeltaP").draw(gr31);

        H2F h32 = dg_dp.getH2F("hi-deltap-phi2");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f32 = new F1D("f32","[amp]*gaus(x,[mean],[sigma])", -0.150, 0.150);
        GraphErrors gr32 = new GraphErrors();
        sliceFit(h32, f32, gr32);
        //DataFitter.fit(fn_deltap_phi2, gr32, "Q");
        myCanvas.getCanvas("DeltaP").cd(10);
        myCanvas.getCanvas("DeltaP").draw(gr32);

        H2F h33 = dg_dp.getH2F("hi-deltap-phi3");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f33 = new F1D("f33","[amp]*gaus(x,[mean],[sigma])", -0.150, 0.150);
        GraphErrors gr33 = new GraphErrors();
        sliceFit(h33, f33, gr33);
        //DataFitter.fit(fn_deltap_phi3, gr33, "Q");
        myCanvas.getCanvas("DeltaP").cd(11);
        myCanvas.getCanvas("DeltaP").draw(gr33);

        H2F h4 = dg_dth.getH2F("hi-deltatheta-p");
        //H2F h1 = hi_fkpckp_deltap_p;
        F1D f4 = new F1D("f4","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        GraphErrors gr4 = new GraphErrors();
        sliceFit(h4, f4, gr4);
        //DataFitter.fit(fn_deltatheta_p, gr4, "Q");
        myCanvas.getCanvas("DeltaTheta").cd(6);
        myCanvas.getCanvas("DeltaTheta").draw(gr4);

        H2F h5 = dg_dth.getH2F("hi-deltatheta-theta");
        //H2F h2 = hi_fkpckp_deltap_theta;
        F1D f5 = new F1D("f5","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        GraphErrors gr5 = new GraphErrors();
        sliceFit(h5, f5, gr5);
        //DataFitter.fit(fn_deltatheta_theta, gr5, "Q");
        myCanvas.getCanvas("DeltaTheta").cd(7);
        myCanvas.getCanvas("DeltaTheta").draw(gr5);

        H2F h6 = dg_dth.getH2F("hi-deltatheta-phi");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f6 = new F1D("f6","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        GraphErrors gr6 = new GraphErrors();
        sliceFit(h6, f6, gr6);
        myCanvas.getCanvas("DeltaTheta").cd(8);
        myCanvas.getCanvas("DeltaTheta").draw(gr6);

        H2F h61 = dg_dth.getH2F("hi-deltatheta-phi1");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f61 = new F1D("f61","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        GraphErrors gr61 = new GraphErrors();
        sliceFit(h61, f61, gr61);
        //DataFitter.fit(fn_deltatheta_phi1, gr61, "Q");
        myCanvas.getCanvas("DeltaTheta").cd(9);
        myCanvas.getCanvas("DeltaTheta").draw(gr61);

        H2F h62 = dg_dth.getH2F("hi-deltatheta-phi2");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f62 = new F1D("f62","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        GraphErrors gr62 = new GraphErrors();
        sliceFit(h62, f62, gr62);
        //DataFitter.fit(fn_deltatheta_phi2, gr62, "Q");
        myCanvas.getCanvas("DeltaTheta").cd(10);
        myCanvas.getCanvas("DeltaTheta").draw(gr62);

        H2F h63 = dg_dth.getH2F("hi-deltatheta-phi3");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f63 = new F1D("f63","[amp]*gaus(x,[mean],[sigma])", -10.0, 10.0);
        GraphErrors gr63 = new GraphErrors();
        sliceFit(h63, f63, gr63);
        //DataFitter.fit(fn_deltatheta_phi3, gr63, "Q");
        myCanvas.getCanvas("DeltaTheta").cd(11);
        myCanvas.getCanvas("DeltaTheta").draw(gr63);

        H2F h7 = dg_dphi.getH2F("hi-deltaphi-p");
        //H2F h1 = hi_fkpckp_deltap_p;
        F1D f7 = new F1D("f7","[amp]*gaus(x,[mean],[sigma])", -4.0, 4.0);
        GraphErrors gr7 = new GraphErrors();
        sliceFit(h7, f7, gr7);
        myCanvas.getCanvas("DeltaPhi").cd(6);
        myCanvas.getCanvas("DeltaPhi").draw(gr7);

        H2F h8 = dg_dphi.getH2F("hi-deltaphi-theta");
        //H2F h2 = hi_fkpckp_deltap_theta;
        F1D f8 = new F1D("f8","[amp]*gaus(x,[mean],[sigma])", -4.0, 4.0);
        GraphErrors gr8 = new GraphErrors();
        sliceFit(h8, f8, gr8);
        myCanvas.getCanvas("DeltaPhi").cd(7);
        myCanvas.getCanvas("DeltaPhi").draw(gr8);

        H2F h9 = dg_dphi.getH2F("hi-deltaphi-phi");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f9 = new F1D("f9","[amp]*gaus(x,[mean],[sigma])", -4.0, 4.0);
        GraphErrors gr9 = new GraphErrors();
        sliceFit(h9, f9, gr9);
        myCanvas.getCanvas("DeltaPhi").cd(8);
        myCanvas.getCanvas("DeltaPhi").draw(gr9);

        H2F h91 = dg_dphi.getH2F("hi-deltaphi-phi1");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f91 = new F1D("f91","[amp]*gaus(x,[mean],[sigma])", -4.0, 4.0);
        GraphErrors gr91 = new GraphErrors();
        sliceFit(h91, f91, gr91);
        //DataFitter.fit(fn_deltaphi_phi1, gr91, "Q");
        myCanvas.getCanvas("DeltaPhi").cd(9);
        myCanvas.getCanvas("DeltaPhi").draw(gr91);

        H2F h92 = dg_dphi.getH2F("hi-deltaphi-phi2");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f92 = new F1D("f92","[amp]*gaus(x,[mean],[sigma])", -4.0, 4.0);
        GraphErrors gr92 = new GraphErrors();
        sliceFit(h92, f92, gr92);
        //DataFitter.fit(fn_deltaphi_phi2, gr92, "Q");
        myCanvas.getCanvas("DeltaPhi").cd(10);
        myCanvas.getCanvas("DeltaPhi").draw(gr92);

        H2F h93 = dg_dphi.getH2F("hi-deltaphi-phi3");
        //H2F h3 = hi_fkpckp_deltap_phi
        F1D f93 = new F1D("f93","[amp]*gaus(x,[mean],[sigma])", -4.0, 4.0);
        GraphErrors gr93 = new GraphErrors();
        sliceFit(h93, f93, gr93);
        //DataFitter.fit(fn_deltaphi_phi3, gr93, "Q");
        myCanvas.getCanvas("DeltaPhi").cd(11);
        myCanvas.getCanvas("DeltaPhi").draw(gr93);

    }

    public void plot(){
        myCanvas = new EmbeddedCanvasTabbed("Kinematics","All","CPIP-missing", "DeltaP","DeltaTheta" ,"DeltaPhi");


        myCanvas.getCanvas("Kinematics").divide(2, 2);
        myCanvas.getCanvas("Kinematics").setSize(1600, 1000);
        myCanvas.getCanvas("Kinematics").setGridX(false);
        myCanvas.getCanvas("Kinematics").setGridY(false);
        myCanvas.getCanvas("Kinematics").setAxisFontSize(18);
        myCanvas.getCanvas("Kinematics").setAxisTitleSize(24);
        myCanvas.getCanvas("Kinematics").draw(dg_kinematics);
        myCanvas.getCanvas("Kinematics").getPad(2).getAxisZ().setLog(true);

        myCanvas.getCanvas("All").divide(2, 2);
        myCanvas.getCanvas("All").setSize(1600, 1000);
        myCanvas.getCanvas("All").setGridX(false);
        myCanvas.getCanvas("All").setGridY(false);
        myCanvas.getCanvas("All").setAxisFontSize(18);
        myCanvas.getCanvas("All").setAxisTitleSize(24);
        myCanvas.getCanvas("All").draw(dg_all);

        myCanvas.getCanvas("CPIP-missing").divide(3, 2);
        myCanvas.getCanvas("CPIP-missing").setSize(1600, 1000);
        myCanvas.getCanvas("CPIP-missing").setGridX(false);
        myCanvas.getCanvas("CPIP-missing").setGridY(false);
        myCanvas.getCanvas("CPIP-missing").setAxisFontSize(18);
        myCanvas.getCanvas("CPIP-missing").setAxisTitleSize(24);
        myCanvas.getCanvas("CPIP-missing").draw(dg_dpdthdphi);

        myCanvas.getCanvas("DeltaP").divide(6, 2);
        myCanvas.getCanvas("DeltaP").setSize(1600, 1000);
        myCanvas.getCanvas("DeltaP").setGridX(false);
        myCanvas.getCanvas("DeltaP").setGridY(false);
        myCanvas.getCanvas("DeltaP").setAxisFontSize(18);
        myCanvas.getCanvas("DeltaP").setAxisTitleSize(24);
        myCanvas.getCanvas("DeltaP").draw(dg_dp);

        myCanvas.getCanvas("DeltaP").getPad(6).getAxisY().setRange(-0.1, 0.1);//-0.05, 0.05
        myCanvas.getCanvas("DeltaP").getPad(7).getAxisY().setRange(-0.1, 0.1);
        myCanvas.getCanvas("DeltaP").getPad(8).getAxisY().setRange(-0.1, 0.1);
        myCanvas.getCanvas("DeltaP").getPad(9).getAxisY().setRange(-0.1, 0.1);
        myCanvas.getCanvas("DeltaP").getPad(10).getAxisY().setRange(-0.1, 0.1);
        myCanvas.getCanvas("DeltaP").getPad(11).getAxisY().setRange(-0.1, 0.1);
        

        myCanvas.getCanvas("DeltaTheta").divide(6, 2);
        myCanvas.getCanvas("DeltaTheta").setSize(1600, 1000);
        myCanvas.getCanvas("DeltaTheta").setGridX(false);
        myCanvas.getCanvas("DeltaTheta").setGridY(false);
        myCanvas.getCanvas("DeltaTheta").setAxisFontSize(18);
        myCanvas.getCanvas("DeltaTheta").setAxisTitleSize(24);
        myCanvas.getCanvas("DeltaTheta").draw(dg_dth);
        /*
        myCanvas.getCanvas("DeltaTheta").getPad(6).getAxisY().setRange(-1.50, 1.50);
        myCanvas.getCanvas("DeltaTheta").getPad(7).getAxisY().setRange(-1.50, 1.50);
        myCanvas.getCanvas("DeltaTheta").getPad(8).getAxisY().setRange(-1.50, 1.50);
        myCanvas.getCanvas("DeltaTheta").getPad(9).getAxisY().setRange(-1.50, 1.50);
        myCanvas.getCanvas("DeltaTheta").getPad(10).getAxisY().setRange(-1.50, 1.50);
        myCanvas.getCanvas("DeltaTheta").getPad(11).getAxisY().setRange(-1.50, 1.50);
        */
     //   /*
        myCanvas.getCanvas("DeltaTheta").getPad(6).getAxisY().setRange(-3.0, 3.0);
        myCanvas.getCanvas("DeltaTheta").getPad(7).getAxisY().setRange(-3.0, 3.0);
        myCanvas.getCanvas("DeltaTheta").getPad(8).getAxisY().setRange(-3.0, 3.0);
        myCanvas.getCanvas("DeltaTheta").getPad(9).getAxisY().setRange(-3.0, 3.0);
        myCanvas.getCanvas("DeltaTheta").getPad(10).getAxisY().setRange(-3.0, 3.0);
        myCanvas.getCanvas("DeltaTheta").getPad(11).getAxisY().setRange(-3.0, 3.0);
       // */

        myCanvas.getCanvas("DeltaPhi").divide(6, 2);
        myCanvas.getCanvas("DeltaPhi").setSize(1600, 1000);
        myCanvas.getCanvas("DeltaPhi").setGridX(false);
        myCanvas.getCanvas("DeltaPhi").setGridY(false);
        myCanvas.getCanvas("DeltaPhi").setAxisFontSize(18);
        myCanvas.getCanvas("DeltaPhi").setAxisTitleSize(24);
        myCanvas.getCanvas("DeltaPhi").draw(dg_dphi);
       // /*
        myCanvas.getCanvas("DeltaPhi").getPad(6).getAxisY().setRange(-1.0, 1.0);
        myCanvas.getCanvas("DeltaPhi").getPad(7).getAxisY().setRange(-1.0, 1.0);
        myCanvas.getCanvas("DeltaPhi").getPad(8).getAxisY().setRange(-1.0, 1.0);
        myCanvas.getCanvas("DeltaPhi").getPad(9).getAxisY().setRange(-1.0, 1.0);
        myCanvas.getCanvas("DeltaPhi").getPad(10).getAxisY().setRange(-1.0, 1.0);
        myCanvas.getCanvas("DeltaPhi").getPad(11).getAxisY().setRange(-1.0, 1.0);
        //*/

    }
    public void saveCanvas(){
        myCanvas.getCanvas("Kinematics").save("eppippim_kinematics.png");
        myCanvas.getCanvas("All").save("eppippim_missing_all.png");
        myCanvas.getCanvas("CPIP-missing").save("eppippim_missing_cpip.png");
        myCanvas.getCanvas("DeltaP").save("cpip_deltaP.png");
        myCanvas.getCanvas("DeltaTheta").save("cpip_deltaTheta.png");
        myCanvas.getCanvas("DeltaPhi").save("cpip_deltaPhi.png");

    }

    public void showplots() {
        JFrame frame = new JFrame("eppippim");
        frame.setSize(1600, 1000);
        frame.add(myCanvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    public static void main(String[] args) {
        //System.setProperty("java.awt.headless", "true");
        GStyle.setPalette("kRainBow");
        GStyle.getH1FAttributes().setOptStat("1110");
        //GStyle.getFunctionAttributes().setOptStat("1100");

        int count = 0;
        //int maxevents = 100000;
        long maxevents = 200000000000;
        eppippim ana = new eppippim();
        //eppippim ana = new eppippim(10.604f);
        System.out.println(String.format(">>> files from list %s >>>", args[0]));
        String filelist;// = "list_of_files.txt";
        filelist = args[0];
        List<String> toProcessFileNames = new ArrayList<String>();
        File file = new File(filelist);
        Scanner read;
        try {
                read = new Scanner(file);
                do {
                        String filename = read.next();
                        toProcessFileNames.add(filename);

                }while (read.hasNext());
                read.close();
        }catch(IOException e){
                e.printStackTrace();
        }
        int progresscount=0;
        int filetot = toProcessFileNames.size();
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        long startTime = System.currentTimeMillis();
        long previousTime = System.currentTimeMillis();
        /*
        HipoDataSync  fileWriter1 = new HipoDataSync();
        HipoDataSync  fileWriter2 = new HipoDataSync();
        fileWriter1.open("filter_ecpipfpimfp.hipo");
        fileWriter2.open("filter_efpipcpimfp.hipo");
        // */



        for (String runstrg : toProcessFileNames) if(count < maxevents){
            progresscount++;
            System.out.println(String.format(">>>> file %s >>>>", runstrg));
            File varTmpDir = new File(runstrg);
            if(!varTmpDir.exists()){System.out.println("FILE DOES NOT EXIST");continue;}
            System.out.println("READING NOW " + runstrg);
            HipoDataSource reader = new HipoDataSource();
            reader.open(runstrg);
            int filecount = 0;
            while(reader.hasEvent() && count < maxevents) {
                //DataEvent event = reader.getNextEvent();
                HipoDataEvent event = reader.getNextEvent();
                ana.processEvent(event);

                /*
                if(ana.select_ecpipfpimfp()) {
                    fileWriter1.writeEvent(event);
                }

                if(ana.select_efpipcpimfp()) {
                    fileWriter2.writeEvent(event);
                }
                //*/
                
                filecount++;count++;
                if(count%1000000 == 0){
                    long nowTime = System.currentTimeMillis();
                    long elapsedTime = nowTime - previousTime;
                    long totalTime = nowTime - startTime;
                    elapsedTime = elapsedTime/1000;
                    totalTime = totalTime/1000;
                    Date date = new Date();
                    String TimeString = "          time : " + dateFormat.format(date) + " , last elapsed : " + elapsedTime + "s ; total elapsed : " + totalTime + "s";
                    String diagnost = String.format(" file %d/%d", progresscount, filetot);
                    System.out.println(diagnost + TimeString); 
                    previousTime = nowTime;
                    
                }               
                
            }          
            reader.close();
        }
        /*
        fileWriter1.close();
        fileWriter2.close();
        //*/        
        System.out.println("Total events : " + count);
        System.out.println("fn_deltap_phi1 Before fitting : " + ana.fn_deltap_phi1.getParameter(0));
        //ana.analyze();
        ana.plot();
        ana.plotGraph();
        ana.showplots();
        //ana.saveCanvas();
        // parameter after fitting function
        //System.out.println("fn_deltap_phi1 After fitting : " + ana.fn_deltap_phi1.getParameter(0));
        //ana.plot();
        //ana.save();
        
        System.out.println("Good Bye !!!!!");

    }



 }

