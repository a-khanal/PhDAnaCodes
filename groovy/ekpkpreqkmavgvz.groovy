import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;
import org.jlab.jnp.physics.*;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.pdg.PDGParticle;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.LatexText;
import org.jlab.groot.ui.LatexTextTools;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import javax.swing.JFrame;

/**
 * @author akhanal
 *
 */

public class ekpkpreqkmavgvz {
	//public int NFTElec;
	
	public float  Eb, Mp;
	public float rec_EB, rec_beam_px, rec_beam_py, rec_beam_pz;
	//public float EB, Eb, Mp;
	public float STT, RFT, FTSTT, vt;
	
	public Particle Vprot, Vpip, Vpim, Vkp, Vkm;
	public Particle Vprotc, Vpipc, Vpimc, Vkpc, Vkmc;
	
	public LorentzVector Ve_consCorn;
	public LorentzVector VB, VT, Ve, VGS, VhadronSystm, Vpim_correct, Vkp_correct;

	public Vector3D e_ftCal_hitPosition;
	
	public boolean found_eFT, found_Lambda, found_Sigma, found_Cascade, found_recPip;
	public int e_ft_part_ind;
	
	public float e_mom, e_the, e_phi;
	public float e_xB, e_Q2, e_W, e_virphoton;
	
	public int prot_part_ind, prot_FTOF_pad1b;
	public float prot_mom, prot_the, prot_phi, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz, prot_ftb_beta, prot_FTOF1b_t, prot_FTOF1b_path, prot_FTOF1b_vt;
	public int prot_CTOF_pad;
	public float prot_CTOF_t, prot_CTOF_path, prot_CTOF_vt;
	
	public int pip_part_ind, pip_FTOF_pad1b;
	public float pip_mom, pip_the, pip_phi, pip_px, pip_py, pip_pz,  pip_vx, pip_vy, pip_vz, pip_status, pip_ftb_beta, pip_FTOF1b_t, pip_FTOF1b_path, pip_FTOF1b_vt;
	public int pip_CTOF_pad;
	public float pip_CTOF_t, pip_CTOF_path, pip_CTOF_vt;
	
	
	public int pim_part_ind, pim_FTOF_pad1b;
	public float pim_mom, pim_the, pim_phi, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz, pim_ftb_beta, pim_FTOF1b_t, pim_FTOF1b_path, pim_FTOF1b_vt;
	public int pim_CTOF_pad;
	public float pim_CTOF_t, pim_CTOF_path, pim_CTOF_vt;
	
	public int kp_part_ind, kp_FTOF_pad1b;
	public float kp_mom, kp_the, kp_phi, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz, kp_ftb_beta, kp_FTOF1b_t, kp_FTOF1b_path, kp_FTOF1b_vt;
	public int kp_CTOF_pad;
	public float kp_CTOF_t, kp_CTOF_path, kp_CTOF_vt;
	
	public int km_part_ind, km_FTOF_pad1b;
	public float km_mom, km_the, km_phi, km_px, km_py, km_pz, km_vx, km_vy, km_vz, km_ftb_beta, km_FTOF1b_t, km_FTOF1b_path, km_FTOF1b_vt;
	public int km_CTOF_pad;
	public float km_CTOF_t, km_CTOF_path, km_CTOF_vt;
	
	
	
	public float ekpkp_MM_req_pim, efkpckp_MM_req_pim, efkpckp_MM_req_cpim, efkpckp_MM_req_km, efkpckp_MM_req_ckm;
	public F1D F_prot_beta_mom, F_kp_beta_mom, F_pip_beta_mom;
	
	
	public DataLine L_lamda, L_sigma;
	
	
	public H1F H_FT_W, H_FT_Q2, H_virphoton;
	public H1F  H_ekpkp_MM_req_pim, H_efkpckp_MM_req_pim, H_efkpckp_MM_req_cpim, H_efkpckp_MM_req_km, H_efkpckp_MM_req_ckm;
	// particles vz, path_FTOF1b, time_FTOF1b
	
	public H1F H_kp_vz, H_km_vz, H_prot_vz, H_pip_vz, H_pim_vz;
	public H1F H_pipc_vz, H_pimc_vz, H_kpc_vz, H_kmc_vz, H_protc_vz; 
	public H1F H_prot_FTOF1b_path, H_pip_FTOF1b_path,  H_pim_FTOF1b_path, H_kp_FTOF1b_path, H_km_FTOF1b_path;
	public H1F H_prot_FTOF1b_t, H_pip_FTOF1b_t,  H_pim_FTOF1b_t, H_kp_FTOF1b_t, H_km_FTOF1b_t;
	public H1F H_pip_CTOF_path, H_pim_CTOF_path, H_prot_CTOF_path,  H_kp_CTOF_path, H_km_CTOF_path;
	public H1F H_pim_CTOF_t, H_pip_CTOF_t, H_prot_CTOF_t, H_kp_CTOF_t, H_km_CTOF_t;
	

	public H2F H_pip_vt_p, H_prot_vt_p, H_pim_vt_p, H_kp_vt_p, H_km_vt_p;
	public H2F H_pipc_vt_p, H_pimc_vt_p, H_kpc_vt_p, H_kmc_vt_p, H_protc_vt_p;
	
	public H1F hi_fd_pos_mass, hi_fd_neg_mass, hi_cd_pos_mass, hi_cd_neg_mass;
	public H2F H_FD_pos_beta_mom, H_FD_neg_beta_mom, H_FD_neutral_beta_mom;
	public H2F H_FD_pos_mass_mom, H_FD_neg_mass_mom, H_FD_neutral_mass_mom;
	public H2F H_FD_pos_mass_the, H_FD_neg_mass_the, H_FD_neutral_mass_the;
	public H2F H_FD_pos_mass_phi, H_FD_neg_mass_phi, H_FD_neutral_mass_phi;
	public H2F H_CD_pos_beta_mom, H_CD_neg_beta_mom, H_CD_neutral_beta_mom;
	public H2F H_CD_pos_mass_mom, H_CD_neg_mass_mom, H_CD_neutral_mass_mom;
	public H2F H_CD_pos_mass_the, H_CD_neg_mass_the, H_CD_neutral_mass_the;
	public H2F H_CD_pos_mass_phi, H_CD_neg_mass_phi, H_CD_neutral_mass_phi;

	public List<Particle> fpips, fpims, fkps, fkms, fprots; // = new ArrayList<Particle>();
	public List<Particle> cpips, cpims, ckps, ckms, cprots;
	public List<Particle> ftels, pips, pims, kps, kms, prots;
	public List<Particle> mcels, mckps, mckms;

    public Particle mckp1;
    public Particle mckp2;
	// for e kp kp km detected
	public Particle reckp1;
    public Particle reckp2;
	public float rec_kp1_p, rec_kp1_the, rec_kp1_phi, rec_kp1_vz, rec_kp2_p, rec_kp2_the, rec_kp2_phi, rec_kp2_vz;
	public float rec_km_p, rec_km_the, rec_km_phi, rec_km_vz;
	public float ekpkpkm_MM_ekpkp, ekpkpkm_MM_ekpkpkm, ekpkpkm_IM_kmlambda, ekpkpkm_IM_kmsigma;
	public float ekpkpkm_MM_ekpkp_withcorrn, ekpkpkm_MM_ekpkpkm_withcorrn;
	public float ekpkpkmprot_MM2;

	
	public F1D fn_rec_e_dp, fn_rec_e_dtheta, fn_rec_e_dphi, fn_rec_e_dvx, fn_rec_e_dvy, fn_rec_e_dvz;
	public F1D fn_rec_kp1_dp, fn_rec_kp1_dtheta, fn_rec_kp1_dphi, fn_rec_kp1_dvx, fn_rec_kp1_dvy, fn_rec_kp1_dvz;
	public F1D fn_rec_kp2_dp, fn_rec_kp2_dtheta, fn_rec_kp2_dphi, fn_rec_kp2_dvx, fn_rec_kp2_dvy, fn_rec_kp2_dvz;
	public F1D fn_rec_km_dp, fn_rec_km_dtheta, fn_rec_km_dphi, fn_rec_km_dvx, fn_rec_km_dvy, fn_rec_km_dvz;

	public H2F H_FT_e_beta_mom;
	public H2F H_FT_e_t_f, H_FT_e_p_f, H_FT_e_p_the;
	public H1F hi_rec_e_dp, hi_rec_e_dtheta, hi_rec_e_dphi;
	public H2F hi_rec_e_dp_p, hi_rec_e_dp_theta, hi_rec_e_dp_phi, hi_rec_e_dp_vz, hi_rec_e_dtheta_p, hi_rec_e_dtheta_theta, hi_rec_e_dtheta_phi, hi_rec_e_dtheta_vz, hi_rec_e_dphi_p, hi_rec_e_dphi_theta, hi_rec_e_dphi_phi, hi_rec_e_dphi_vz, hi_rec_e_dvz_p, hi_rec_e_dvz_theta, hi_rec_e_dvz_phi, hi_rec_e_dvz_vz;
	public H2F H_FT_W_Q2, H_FT_e_xB_Q2;

	// momentum, theta and phi resolution for kp1, kp2 and km
	public H1F hi_rec_kp1_dp, hi_rec_kp1_dtheta, hi_rec_kp1_dphi;
	public H2F hi_rec_kp1_dp_p, hi_rec_kp1_dp_theta, hi_rec_kp1_dp_phi, hi_rec_kp1_dp_vz, hi_rec_kp1_dtheta_p, hi_rec_kp1_dtheta_theta, hi_rec_kp1_dtheta_phi, hi_rec_kp1_dtheta_vz, hi_rec_kp1_dphi_p, hi_rec_kp1_dphi_theta, hi_rec_kp1_dphi_phi, hi_rec_kp1_dphi_vz, hi_rec_kp1_dvz_p, hi_rec_kp1_dvz_theta, hi_rec_kp1_dvz_phi, hi_rec_kp1_dvz_vz;
	public H1F hi_rec_kp2_dp, hi_rec_kp2_dtheta, hi_rec_kp2_dphi;
	public H2F hi_rec_kp2_dp_p, hi_rec_kp2_dp_theta, hi_rec_kp2_dp_phi, hi_rec_kp2_dp_vz, hi_rec_kp2_dtheta_p, hi_rec_kp2_dtheta_theta, hi_rec_kp2_dtheta_phi, hi_rec_kp2_dtheta_vz, hi_rec_kp2_dphi_p, hi_rec_kp2_dphi_theta, hi_rec_kp2_dphi_phi, hi_rec_kp2_dphi_vz, hi_rec_kp2_dvz_p, hi_rec_kp2_dvz_theta, hi_rec_kp2_dvz_phi, hi_rec_kp2_dvz_vz;
	public H1F hi_rec_km_dp, hi_rec_km_dtheta, hi_rec_km_dphi;
	public H2F hi_rec_km_dp_p, hi_rec_km_dp_theta, hi_rec_km_dp_phi, hi_rec_km_dp_vz, hi_rec_km_dtheta_p, hi_rec_km_dtheta_theta, hi_rec_km_dtheta_phi, hi_rec_km_dtheta_vz, hi_rec_km_dphi_p, hi_rec_km_dphi_theta, hi_rec_km_dphi_phi, hi_rec_km_dphi_vz, hi_rec_km_dvz_p, hi_rec_km_dvz_theta, hi_rec_km_dvz_phi, hi_rec_km_dvz_vz;

	// scatter plotes for hadrons 
	public H2F hi_rec_kp1_p_the, hi_rec_kp2_p_the, hi_rec_km_p_the, hi_mc_e_p_the, hi_mc_kp1_p_the, hi_mc_kp2_p_the, hi_mc_km_p_the;
	public H2F hi_rec_kp1_p_phi, hi_rec_kp2_p_phi, hi_rec_km_p_phi, hi_mc_e_p_phi, hi_mc_kp1_p_phi, hi_mc_kp2_p_phi, hi_mc_km_p_phi;
	public H2F hi_rec_kp1_the_phi, hi_rec_kp2_the_phi, hi_rec_km_the_phi, hi_mc_e_the_phi, hi_mc_kp1_the_phi, hi_mc_kp2_the_phi, hi_mc_km_the_phi;

	//vertex resolution for hadrons 
	public H1F hi_rec_e_dvx, hi_rec_kp1_dvx, hi_rec_kp2_dvx, hi_rec_km_dvx;
	public H1F hi_rec_e_dvy, hi_rec_kp1_dvy, hi_rec_kp2_dvy, hi_rec_km_dvy;
	public H1F hi_rec_e_dvz, hi_rec_kp1_dvz, hi_rec_kp2_dvz, hi_rec_km_dvz;

	//mass spectrum
	public H1F hi_ekpkpkm_MM_ekpkp, hi_ekpkpkm_MM_ekpkpkm, hi_ekpkpkm_IM_kmlambda, hi_ekpkpkm_IM_kmsigma;
	public H1F hi_ekpkpkm_MM_ekpkp_withcorrn, hi_ekpkpkm_MM_ekpkpkm_withcorrn
	public F1D f1_xi, fn_xi_with_cor, f1_mc_xi, f1_lambda, f1_sigma, f_gaus, fn_im_kmlambda, fn_im_kmsigma;
	public H1F hi_mc_ekpkpkm_mm_ekpkp, hi_mc_ekpkpkm_mm_ekpkpkm;


	public H1F hi_ekpkpkmprot_MM2;

	public H1F hi_pip_counter, hi_pim_counter, hi_kp_counter, hi_km_counter, hi_prot_counter, hi_fpip_counter, hi_fpim_counter, hi_fkp_counter, hi_fkm_counter, hi_fprot_counter, hi_cpip_counter, hi_cpim_counter, hi_ckp_counter, hi_ckm_counter, hi_cprot_counter;

	public EmbeddedCanvasTabbed myCanvas;
	public DataGroup dg_rec_electron, dg_rec_kp1, dg_rec_kp2, dg_rec_km, dg_rec_xi, dg_rec_p, dg_rec_pim, dg_vtime, dg_fdtof, dg_cdtof, dg_vz, dg_tof_t, dg_tof_path, dg_counter;
	public DataGroup dg_rec_e_resolution, dg_rec_kp1_resolution, dg_rec_kp2_resolution, dg_rec_km_resolution;
	public DataGroup dg_mm;
	public DataGroup dg_req;

	
	
	public ekpkpreqkmavgvz() {

		final int BLUE = 9;
		final int LIGHTGREEN = 3;
		final int LIGHTBROWN = 45;
		final int PINK = 46;
	//	NFTElec = 0;
	//	Eb = 10.575f;
		Eb = 10.604f;
	//	Eb = 7.54626f;
		Mp = (float) PDGDatabase.getParticleById(2212).mass();
	//	Mp = 0.93827f;
		
		VB = new LorentzVector(0, 0, Eb, Eb);
		VT = new LorentzVector(0, 0, 0, Mp);

		
		// theoretical 1D functions for proton, kaon and pion
		
		F_prot_beta_mom = new F1D("F_prot_beta_mom", "x/sqrt(0.93827*0.93827+x*x)", 0.3, 4.0);
		F_prot_beta_mom.setLineWidth(2);
		F_prot_beta_mom.setLineColor(BLUE);
		F_kp_beta_mom = new F1D("F_kp_beta_mom", "x/sqrt(0.49367*0.49367+x*x)", 0.3, 4.0);
		F_kp_beta_mom.setLineColor(BLUE);
		F_kp_beta_mom.setLineWidth(2);
		F_pip_beta_mom = new F1D("F_pip_beta_mom", "x/sqrt(0.13957*0.13957+x*x)", 0.3, 4.0);
		F_pip_beta_mom.setLineColor(BLUE);
		F_pip_beta_mom.setLineWidth(2);
		


	// FT electron overview
		// generated electrons
		hi_mc_e_p_the = new H2F("hi_mc_e_p_the", "hi_mc_e_p_the", 100, 2, 5, 100, 0, 5);
		hi_mc_e_p_the.setTitleX("#theta (^o)");
		hi_mc_e_p_the.setTitleY("p (GeV)");

		hi_mc_e_p_phi = new H2F("hi_mc_e_p_phi", "hi_mc_e_p_phi", 100, -180, 180, 100, 0, 5);
		hi_mc_e_p_phi.setTitleX("#phi (^o)");
		hi_mc_e_p_phi.setTitleY("p (GeV)");

		hi_mc_e_the_phi = new H2F("hi_mc_e_the_phi", "hi_mc_e_the_phi", 100, -180, 180, 100, 2, 5);
		hi_mc_e_the_phi.setTitleX("#phi (^o)");
		hi_mc_e_the_phi.setTitleY("#theta (^o)");

		//reconstructed
		H_FT_e_t_f = new H2F("H_FT_e_t_f", "H_FT_e_t_f", 100, -180, 180, 100, 2, 5);
		H_FT_e_t_f.setTitle("electron #theta vs #phi");
		H_FT_e_t_f.setTitleX("#phi (^o)");
		H_FT_e_t_f.setTitleY("#theta (^o)");
		
		H_FT_e_p_the = new H2F("H_FT_e_p_the", "H_FT_e_p_the", 100, 2, 5, 100, 0, 5);
		H_FT_e_p_the.setTitle("electron p vs #theta (^o)");
		H_FT_e_p_the.setTitleX("#theta (^o)");
		H_FT_e_p_the.setTitleY("p (GeV)");

		H_FT_e_p_f = new H2F("H_FT_e_p_f", "H_FT_e_p_f", 100, -180, 180, 100, 0, 5);
		H_FT_e_p_f.setTitle("electron p vs #phi");
		H_FT_e_p_f.setTitleX("#phi (^o)");
		H_FT_e_p_f.setTitleY("p (GeV)");

		//hi_rec_e_dvx = new H1F("hi_rec_e_dvx", "hi_rec_e_dvx", 100, -10.0, 10.0);
		//hi_rec_e_dvx.setFillColor(LIGHTGREEN);
		//hi_rec_e_dvx.setTitleX("#Delta vx");
		//hi_rec_e_dvy = new H1F("hi_rec_e_dvy", "hi_rec_e_dvy", 100, -10.0, 10.0);
		//hi_rec_e_dvy.setFillColor(LIGHTGREEN);
		//hi_rec_e_dvy.setTitleX("#Delta vy");
		hi_rec_e_dvz = new H1F("hi_rec_e_dvz", "hi_rec_e_dvz", 100, -3.0, 3.0);
		hi_rec_e_dvz.setFillColor(LIGHTGREEN);
		hi_rec_e_dvz.setTitleX("#Delta vz");
		hi_rec_e_dp = new H1F("hi_rec_e_dp", "hi_rec_e_dp", 100, -0.2, 0.2);
		hi_rec_e_dp.setTitleX("#DeltaP/P");
		hi_rec_e_dp.setFillColor(LIGHTGREEN);
		hi_rec_e_dtheta = new H1F("hi_rec_e_dtheta", "hi_rec_e_dtheta", 100, -0.5, 0.5);
		hi_rec_e_dtheta.setTitleX("#Delta#theta (^o)");
		hi_rec_e_dtheta.setFillColor(LIGHTGREEN);
		hi_rec_e_dphi = new H1F("hi_rec_e_dphi","hi_rec_e_dphi", 400, -8.0, 8.0);
		hi_rec_e_dphi.setTitleX("#Delta#phi (^o)");
		hi_rec_e_dphi.setTitleY("count");
		hi_rec_e_dphi.setFillColor(LIGHTGREEN);

		// delta functions 
		fn_rec_e_dp = new F1D("fn_rec_e_dp", "[amp]*gaus(x,[mean],[sigma])", -1, 1);
		fn_rec_e_dp.setLineWidth(2);
		fn_rec_e_dp.setLineColor(2);
		fn_rec_e_dp.setOptStat("1111");
		fn_rec_e_dtheta = new F1D("fn_rec_e_dtheta", "[amp]*gaus(x,[mean],[sigma])", -0.4, 0.4);
		fn_rec_e_dtheta.setLineWidth(2);
		fn_rec_e_dtheta.setLineColor(2);
		fn_rec_e_dtheta.setOptStat("1111");
		fn_rec_e_dphi = new F1D("fn_rec_e_dphi", "[amp]*gaus(x,[mean],[sigma])", -4, 4);
		fn_rec_e_dphi.setLineWidth(2);
		fn_rec_e_dphi.setLineColor(2);
		fn_rec_e_dphi.setOptStat("1111");
		fn_rec_e_dvx = new F1D("fn_rec_e_dvx", "[amp]*gaus(x,[mean],[sigma])", -2.5, 2.5);
		fn_rec_e_dvx.setLineWidth(2);
		fn_rec_e_dvx.setLineColor(2);
		fn_rec_e_dvx.setOptStat("1111");
		fn_rec_e_dvy = new F1D("fn_rec_e_dvy", "[amp]*gaus(x,[mean],[sigma])", -2.5, 2.5);
		fn_rec_e_dvy.setLineWidth(2);
		fn_rec_e_dvy.setLineColor(2);
		fn_rec_e_dvy.setOptStat("1111");
		fn_rec_e_dvz = new F1D("fn_rec_e_dvz", "[amp]*gaus(x,[mean],[sigma])", -2.5, 2.5);
		fn_rec_e_dvz.setLineWidth(2);
		fn_rec_e_dvz.setLineColor(2);
		fn_rec_e_dvz.setOptStat("1111");

		// electron resolution
		hi_rec_e_dp_p =new H2F("hi_rec_e_dp_p","hi_rec_e_dp_p", 100, 0, 5, 100, -0.2, 0.2);
		hi_rec_e_dp_p.setTitleX("p (GeV)");
		hi_rec_e_dp_p.setTitleY("#Delta p/p");
		hi_rec_e_dp_theta =new H2F("hi_rec_e_dp_theta","hi_rec_e_dp_theta", 100, 2, 5, 100, -0.2, 0.2);
		hi_rec_e_dp_theta.setTitleX("#theta (^o)");
		hi_rec_e_dp_theta.setTitleY("#Delta p/p");
		hi_rec_e_dp_phi =new H2F("hi_rec_e_dp_phi","hi_rec_e_dp_phi", 100, 180, -180, 100, -0.2, 0.2);
		hi_rec_e_dp_phi.setTitleX("#phi (^o)");
		hi_rec_e_dp_phi.setTitleY("#Delta p/p");
		hi_rec_e_dp_vz =new H2F("hi_rec_e_dp_vz","hi_rec_e_dp_vz", 100, -10, 10, 100, -0.2, 0.2);
		hi_rec_e_dp_vz.setTitleX("vz");
		hi_rec_e_dp_vz.setTitleY("#Delta p/p");
		hi_rec_e_dtheta_p =new H2F("hi_rec_e_dtheta_p","hi_rec_e_dtheta_p", 100, 0, 5, 100, -0.5, 0.5);
		hi_rec_e_dtheta_p.setTitleX("p (GeV)");
		hi_rec_e_dtheta_p.setTitleY("#Delta#theta (^o)");
		hi_rec_e_dtheta_theta =new H2F("hi_rec_e_dtheta_theta","hi_rec_e_dtheta_theta", 100, 2, 5, 100, -0.5, 0.5);
		hi_rec_e_dtheta_theta.setTitleX("#theta (^o)");
		hi_rec_e_dtheta_theta.setTitleY("#Delta#theta (^o)");
		hi_rec_e_dtheta_phi =new H2F("hi_rec_e_dtheta_phi","hi_rec_e_dtheta_phi", 100, 180, -180, 100, -0.5, 0.5);
		hi_rec_e_dtheta_phi.setTitleX("#phi (^o)");
		hi_rec_e_dtheta_phi.setTitleY("#Delta#theta (^o)");
		hi_rec_e_dtheta_vz =new H2F("hi_rec_e_dtheta_vz","hi_rec_e_dtheta_vz", 100, -10, 10, 100, -0.5, 0.5);
		hi_rec_e_dtheta_vz.setTitleX("vz");
		hi_rec_e_dtheta_vz.setTitleY("#Delta#theta (^o)");
		hi_rec_e_dphi_p =new H2F("hi_rec_e_dphi_p","hi_rec_e_dphi_p", 100, 0, 5, 400, -8.0, 8.0);
		hi_rec_e_dphi_p.setTitleX("p (GeV)");
		hi_rec_e_dphi_p.setTitleY("#Delta#phi (^o)");
		hi_rec_e_dphi_theta =new H2F("hi_rec_e_dphi_theta","hi_rec_e_dphi_theta", 100, 2, 5, 400, -8.0, 8.0);
		hi_rec_e_dphi_theta.setTitleX("#theta (^o)");
		hi_rec_e_dphi_theta.setTitleY("#Delta#phi (^o)");
		hi_rec_e_dphi_phi =new H2F("hi_rec_e_dphi_phi","hi_rec_e_dphi_phi", 100, 180, -180, 400, -8.0, 8.0);
		hi_rec_e_dphi_phi.setTitleX("#phi (^o)");
		hi_rec_e_dphi_phi.setTitleY("#Delta#phi (^o)");
		hi_rec_e_dphi_vz =new H2F("hi_rec_e_dphi_vz","hi_rec_e_dphi_vz", 100, -10, 10, 400, -8.0, 8.0);
		hi_rec_e_dphi_vz.setTitleX("vz");
		hi_rec_e_dphi_vz.setTitleY("#Delta#phi (^o)");
		hi_rec_e_dvz_p =new H2F("hi_rec_e_dvz_p","hi_rec_e_dvz_p", 100, 0, 5, 100, -3.0, 3.0);
		hi_rec_e_dvz_p.setTitleX("p (GeV)");
		hi_rec_e_dvz_p.setTitleY("#Delta vz");
		hi_rec_e_dvz_theta =new H2F("hi_rec_e_dvz_theta","hi_rec_e_dvz_theta", 100, 2, 5, 100, -3.0, 3.0);
		hi_rec_e_dvz_theta.setTitleX("#theta (^o)");
		hi_rec_e_dvz_theta.setTitleY("#Delta vz");
		hi_rec_e_dvz_phi =new H2F("hi_rec_e_dvz_phi","hi_rec_e_dvz_phi", 100, 180, -180, 100, -3.0, 3.0);
		hi_rec_e_dvz_phi.setTitleX("#phi (^o)");
		hi_rec_e_dvz_phi.setTitleY("#Delta vz");
		hi_rec_e_dvz_vz =new H2F("hi_rec_e_dvz_vz","hi_rec_e_dvz_vz", 100, -10, 10, 100, -3.0, 3.0);
		hi_rec_e_dvz_vz.setTitleX("vz");
		hi_rec_e_dvz_vz.setTitleY("#Delta vz");
		
		H_FT_W_Q2 = new H2F("H_FT_W_Q2", "H_FT_W_Q2", 100, 3, 5, 100, 0.00001, 0.5);
		H_FT_W_Q2.setTitle("FT Q^2 vs W");
		H_FT_W_Q2.setTitleX("W ( GeV)");
		H_FT_W_Q2.setTitleY("Q^2 (GeV^2)");
		
		
		H_FT_W = new H1F("H_FT_W", "H_FT_W", 100, 3, 5);
		H_FT_W.setTitle("electron W");
		H_FT_W.setTitleX("W (GeV)");
		H_FT_W.setTitleY("count");
		H_FT_W.setFillColor(LIGHTGREEN);


		H_FT_Q2 = new H1F("H_FT_Q2", "H_FT_Q2", 100, 0.00001, 0.5);
		H_FT_Q2.setFillColor(LIGHTGREEN);
		H_FT_Q2.setTitleX("Q^2 (GeV^2)");

		H_virphoton = new H1F("H_virphoton", "H_virphoton", 100, 0, 12);
		H_virphoton.setFillColor(LIGHTGREEN);
		H_virphoton.setTitleX("E_#gamma (GeV)");

	// KP1
		// generated 
		hi_mc_kp1_p_the = new H2F("hi_mc_kp1_p_the", "hi_mc_kp1_p_the", 100, 0, 100, 100, 0, 8);
		hi_mc_kp1_p_the.setTitleX("#theta (^o)");
		hi_mc_kp1_p_the.setTitleY("p (GeV)");
		hi_mc_kp1_p_phi = new H2F("hi_mc_kp1_p_phi", "hi_mc_kp1_p_phi", 100, -180, 180, 100, 0, 8);
		hi_mc_kp1_p_phi.setTitleX("#phi (^o)");
		hi_mc_kp1_p_phi.setTitleY("p (GeV)");
		hi_mc_kp1_the_phi = new H2F("hi_mc_kp1_the_phi", "hi_mc_kp1_the_phi", 100, -180, 180, 100, 0, 100);
		hi_mc_kp1_the_phi.setTitleX("#phi (^o)");
		hi_mc_kp1_the_phi.setTitleY("#theta (^o)");

		// reconstructed
		hi_rec_kp1_p_the = new H2F("hi_rec_kp1_p_the", "hi_rec_kp1_p_the", 100, 0, 100, 100, 0, 8);
		hi_rec_kp1_p_the.setTitleX("#theta (^o)");
		hi_rec_kp1_p_the.setTitleY("p (GeV)");
		hi_rec_kp1_p_phi = new H2F("hi_rec_kp1_p_phi", "hi_rec_kp1_p_phi", 100, -180, 180, 100, 0, 8);
		hi_rec_kp1_p_phi.setTitleX("#phi (^o)");
		hi_rec_kp1_p_phi.setTitleY("p (GeV)");
		hi_rec_kp1_the_phi = new H2F("hi_rec_kp1_the_phi", "hi_rec_kp1_the_phi", 100, -180, 180, 100, 0, 100);
		hi_rec_kp1_the_phi.setTitleX("#phi (^o)");
		hi_rec_kp1_the_phi.setTitleY("#theta (^o)");


		hi_rec_kp1_dvx = new H1F("hi_rec_kp1_dvx", "hi_rec_kp1_dvx", 100, -3.0, 3.0);
		hi_rec_kp1_dvx.setFillColor(LIGHTGREEN);
		hi_rec_kp1_dvx.setTitleX("#Delta vx");
		hi_rec_kp1_dvy = new H1F("hi_rec_kp1_dvy", "hi_rec_kp1_dvy", 100, -3.0, 3.0);
		hi_rec_kp1_dvy.setFillColor(LIGHTGREEN);
		hi_rec_kp1_dvy.setTitleX("#Delta vy");
		hi_rec_kp1_dvz = new H1F("hi_rec_kp1_dvz", "hi_rec_kp1_dvz", 100, -3.0, 3.0);
		hi_rec_kp1_dvz.setFillColor(LIGHTGREEN);
		hi_rec_kp1_dvz.setTitleX("#Delta vz");
		hi_rec_kp1_dp = new H1F("hi_rec_kp1_dp", "hi_rec_kp1_dp", 100, -0.05, 0.05);
		hi_rec_kp1_dp.setTitleX("#Delta P/P");
		hi_rec_kp1_dp.setFillColor(LIGHTGREEN);
		hi_rec_kp1_dtheta = new H1F("hi_rec_kp1_dtheta", "hi_rec_kp1_dtheta", 100, -0.5, 0.5);
		hi_rec_kp1_dtheta.setTitleX("#Delta#theta (^o)");
		hi_rec_kp1_dtheta.setFillColor(LIGHTGREEN);
		hi_rec_kp1_dphi = new H1F("hi_rec_kp1_dphi","hi_rec_kp1_dphi", 100, -2.0, 2.0);
		hi_rec_kp1_dphi.setTitleX("#Delta#phi (^o)");
		hi_rec_kp1_dphi.setTitleY("count");
		hi_rec_kp1_dphi.setFillColor(LIGHTGREEN);

		// delta functions 
		fn_rec_kp1_dp = new F1D("fn_rec_kp1_dp", "[amp]*gaus(x,[mean],[sigma])", -1, 1);
		fn_rec_kp1_dp.setLineWidth(2);
		fn_rec_kp1_dp.setLineColor(2);
		fn_rec_kp1_dp.setOptStat("1111");
		fn_rec_kp1_dtheta = new F1D("fn_rec_kp1_dtheta", "[amp]*gaus(x,[mean],[sigma])", -0.5, 0.5);
		fn_rec_kp1_dtheta.setLineWidth(2);
		fn_rec_kp1_dtheta.setLineColor(2);
		fn_rec_kp1_dtheta.setOptStat("1111");
		fn_rec_kp1_dphi = new F1D("fn_rec_kp1_dphi", "[amp]*gaus(x,[mean],[sigma])", -4, 4);
		fn_rec_kp1_dphi.setLineWidth(2);
		fn_rec_kp1_dphi.setLineColor(2);
		fn_rec_kp1_dphi.setOptStat("1111");
		fn_rec_kp1_dvx = new F1D("fn_rec_kp1_dvx", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_kp1_dvx.setLineWidth(2);
		fn_rec_kp1_dvx.setLineColor(2);
		fn_rec_kp1_dvx.setOptStat("1111");
		fn_rec_kp1_dvy = new F1D("fn_rec_kp1_dvy", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_kp1_dvy.setLineWidth(2);
		fn_rec_kp1_dvy.setLineColor(2);
		fn_rec_kp1_dvy.setOptStat("1111");
		fn_rec_kp1_dvz = new F1D("fn_rec_kp1_dvz", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_kp1_dvz.setLineWidth(2);
		fn_rec_kp1_dvz.setLineColor(2);
		fn_rec_kp1_dvz.setOptStat("1111");

		hi_rec_kp1_dp_p =new H2F("hi_rec_kp1_dp_p","hi_rec_kp1_dp_p", 100, 0, 8, 100, -0.05, 0.05);
		hi_rec_kp1_dp_p.setTitleX("p (GeV)");
		hi_rec_kp1_dp_p.setTitleY("#Delta p/p");
		hi_rec_kp1_dp_theta =new H2F("hi_rec_kp1_dp_theta","hi_rec_kp1_dp_theta", 100, 0, 100, 100, -0.05, 0.05);
		hi_rec_kp1_dp_theta.setTitleX("#theta (^o)");
		hi_rec_kp1_dp_theta.setTitleY("#Delta p/p");
		hi_rec_kp1_dp_phi =new H2F("hi_rec_kp1_dp_phi","hi_rec_kp1_dp_phi", 100, 180, -180, 100, -0.05, 0.05);
		hi_rec_kp1_dp_phi.setTitleX("#phi (^o)");
		hi_rec_kp1_dp_phi.setTitleY("#Delta p/p");
		hi_rec_kp1_dp_vz =new H2F("hi_rec_kp1_dp_vz","hi_rec_kp1_dp_vz", 100, -10, 10, 100, -0.05, 0.05);
		hi_rec_kp1_dp_vz.setTitleX("vz");
		hi_rec_kp1_dp_vz.setTitleY("#Delta p/p");
		hi_rec_kp1_dtheta_p =new H2F("hi_rec_kp1_dtheta_p","hi_rec_kp1_dtheta_p", 100, 0, 8, 100, -0.5, 0.5);
		hi_rec_kp1_dtheta_p.setTitleX("p (GeV)");
		hi_rec_kp1_dtheta_p.setTitleY("#Delta#theta (^o)");
		hi_rec_kp1_dtheta_theta =new H2F("hi_rec_kp1_dtheta_theta","hi_rec_kp1_dtheta_theta", 100, 0, 100, 100, -0.5, 0.5);
		hi_rec_kp1_dtheta_theta.setTitleX("#theta (^o)");
		hi_rec_kp1_dtheta_theta.setTitleY("#Delta#theta (^o)");
		hi_rec_kp1_dtheta_phi =new H2F("hi_rec_kp1_dtheta_phi","hi_rec_kp1_dtheta_phi", 100, 180, -180, 100, -0.5, 0.5);
		hi_rec_kp1_dtheta_phi.setTitleX("#phi (^o)");
		hi_rec_kp1_dtheta_phi.setTitleY("#Delta#theta (^o)");
		hi_rec_kp1_dtheta_vz =new H2F("hi_rec_kp1_dtheta_vz","hi_rec_kp1_dtheta_vz", 100, -10, 10, 100, -0.5, 0.5);
		hi_rec_kp1_dtheta_vz.setTitleX("vz");
		hi_rec_kp1_dtheta_vz.setTitleY("#Delta#theta (^o)");
		hi_rec_kp1_dphi_p =new H2F("hi_rec_kp1_dphi_p","hi_rec_kp1_dphi_p", 100, 0, 8, 100, -2.0, 2.0);
		hi_rec_kp1_dphi_p.setTitleX("p (GeV)");
		hi_rec_kp1_dphi_p.setTitleY("#Delta#phi (^o)");
		hi_rec_kp1_dphi_theta =new H2F("hi_rec_kp1_dphi_theta","hi_rec_kp1_dphi_theta", 100, 0, 100, 100, -2.0, 2.0);
		hi_rec_kp1_dphi_theta.setTitleX("#theta (^o)");
		hi_rec_kp1_dphi_theta.setTitleY("#Delta#phi (^o)");
		hi_rec_kp1_dphi_phi =new H2F("hi_rec_kp1_dphi_phi","hi_rec_kp1_dphi_phi", 100, 180, -180, 100, -2.0, 2.0);
		hi_rec_kp1_dphi_phi.setTitleX("#phi (^o)");
		hi_rec_kp1_dphi_phi.setTitleY("#Delta#phi (^o)");
		hi_rec_kp1_dphi_vz =new H2F("hi_rec_kp1_dphi_vz","hi_rec_kp1_dphi_vz", 100, -10, 10, 100, -2.0, 2.0);
		hi_rec_kp1_dphi_vz.setTitleX("vz");
		hi_rec_kp1_dphi_vz.setTitleY("#Delta#phi (^o)");
		hi_rec_kp1_dvz_p =new H2F("hi_rec_kp1_dvz_p","hi_rec_kp1_dvz_p", 100, 0, 8, 100, -3.0, 3.0);
		hi_rec_kp1_dvz_p.setTitleX("p (GeV)");
		hi_rec_kp1_dvz_p.setTitleY("#Delta vz");
		hi_rec_kp1_dvz_theta =new H2F("hi_rec_kp1_dvz_theta","hi_rec_kp1_dvz_theta", 100, 0, 100, 100, -3.0, 3.0);
		hi_rec_kp1_dvz_theta.setTitleX("#theta (^o)");
		hi_rec_kp1_dvz_theta.setTitleY("#Delta vz");
		hi_rec_kp1_dvz_phi =new H2F("hi_rec_kp1_dvz_phi","hi_rec_kp1_dvz_phi", 100, 180, -180, 100, -3.0, 3.0);
		hi_rec_kp1_dvz_phi.setTitleX("#phi (^o)");
		hi_rec_kp1_dvz_phi.setTitleY("#Delta vz");
		hi_rec_kp1_dvz_vz =new H2F("hi_rec_kp1_dvz_vz","hi_rec_kp1_dvz_vz", 100, -10, 10, 100, -3.0, 3.0);
		hi_rec_kp1_dvz_vz.setTitleX("vz");
		hi_rec_kp1_dvz_vz.setTitleY("#Delta vz");

	// KP2
		// generated 
		hi_mc_kp2_p_the = new H2F("hi_mc_kp2_p_the", "hi_mc_kp2_p_the", 100, 0, 100, 100, 0, 8);
		hi_mc_kp2_p_the.setTitleX("#theta (^o)");
		hi_mc_kp2_p_the.setTitleY("p (GeV)");
		hi_mc_kp2_p_phi = new H2F("hi_mc_kp2_p_phi", "hi_mc_kp2_p_phi", 100, -180, 180, 100, 0, 8);
		hi_mc_kp2_p_phi.setTitleX("#phi (^o)");
		hi_mc_kp2_p_phi.setTitleY("p (GeV)");
		hi_mc_kp2_the_phi = new H2F("hi_mc_kp2_the_phi", "hi_mc_kp2_the_phi", 100, -180, 180, 100, 0, 100);
		hi_mc_kp2_the_phi.setTitleX("#phi (^o)");
		hi_mc_kp2_the_phi.setTitleY("#theta (^o)");

		// reconstructed
		hi_rec_kp2_p_the = new H2F("hi_rec_kp2_p_the", "hi_rec_kp2_p_the", 100, 0, 100, 100, 0, 8);
		hi_rec_kp2_p_the.setTitleX("#theta (^o)");
		hi_rec_kp2_p_the.setTitleY("p (GeV)");
		hi_rec_kp2_p_phi = new H2F("hi_rec_kp2_p_phi", "hi_rec_kp2_p_phi", 100, -180, 180, 100, 0, 8);
		hi_rec_kp2_p_phi.setTitleX("#phi (^o)");
		hi_rec_kp2_p_phi.setTitleY("p (GeV)");
		hi_rec_kp2_the_phi = new H2F("hi_rec_kp2_the_phi", "hi_rec_kp2_the_phi", 100, -180, 180, 100, 0, 100);
		hi_rec_kp2_the_phi.setTitleX("#phi (^o)");
		hi_rec_kp2_the_phi.setTitleY("#theta (^o)");

		hi_rec_kp2_dvx = new H1F("hi_rec_kp2_dvx", "hi_rec_kp2_dvx", 100, -3.0, 3.0);
		hi_rec_kp2_dvx.setFillColor(LIGHTGREEN);
		hi_rec_kp2_dvx.setTitleX("#Delta vx");
		hi_rec_kp2_dvy = new H1F("hi_rec_kp2_dvy", "hi_rec_kp2_dvy", 100, -3.0, 3.0);
		hi_rec_kp2_dvy.setFillColor(LIGHTGREEN);
		hi_rec_kp2_dvy.setTitleX("#Delta vy");
		hi_rec_kp2_dvz = new H1F("hi_rec_kp2_dvz", "hi_rec_kp2_dvz", 100, -3.0, 3.0);
		hi_rec_kp2_dvz.setFillColor(LIGHTGREEN);
		hi_rec_kp2_dvz.setTitleX("#Delta vz");
		hi_rec_kp2_dp = new H1F("hi_rec_kp2_dp", "hi_rec_kp2_dp", 100, -0.05, 0.05);
		hi_rec_kp2_dp.setTitleX("#delta P/P");
		hi_rec_kp2_dp.setFillColor(LIGHTGREEN);
		hi_rec_kp2_dtheta = new H1F("hi_rec_kp2_dtheta", "hi_rec_kp2_dtheta", 100, -0.5, 0.5);
		hi_rec_kp2_dtheta.setTitleX("#Delta#theta");
		hi_rec_kp2_dtheta.setFillColor(LIGHTGREEN);
		hi_rec_kp2_dphi = new H1F("hi_rec_kp2_dphi","hi_rec_kp2_dphi", 100, -2.0, 2.0);
		hi_rec_kp2_dphi.setTitleX("#Delta#phi");
		hi_rec_kp2_dphi.setTitleY("count");
		hi_rec_kp2_dphi.setFillColor(LIGHTGREEN);

		// delta functions 
		fn_rec_kp2_dp = new F1D("fn_rec_kp2_dp", "[amp]*gaus(x,[mean],[sigma])", -0.1, 0.1);
		fn_rec_kp2_dp.setLineWidth(2);
		fn_rec_kp2_dp.setLineColor(2);
		fn_rec_kp2_dp.setOptStat("1111");
		fn_rec_kp2_dtheta = new F1D("fn_rec_kp2_dtheta", "[amp]*gaus(x,[mean],[sigma])", -0.5, 0.5);
		fn_rec_kp2_dtheta.setLineWidth(2);
		fn_rec_kp2_dtheta.setLineColor(2);
		fn_rec_kp2_dtheta.setOptStat("1111");
		fn_rec_kp2_dphi = new F1D("fn_rec_kp2_dphi", "[amp]*gaus(x,[mean],[sigma])", -4, 4);
		fn_rec_kp2_dphi.setLineWidth(2);
		fn_rec_kp2_dphi.setLineColor(2);
		fn_rec_kp2_dphi.setOptStat("1111");
		fn_rec_kp2_dvx = new F1D("fn_rec_kp2_dvx", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_kp2_dvx.setLineWidth(2);
		fn_rec_kp2_dvx.setLineColor(2);
		fn_rec_kp2_dvx.setOptStat("1111");
		fn_rec_kp2_dvy = new F1D("fn_rec_kp2_dvy", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_kp2_dvy.setLineWidth(2);
		fn_rec_kp2_dvy.setLineColor(2);
		fn_rec_kp2_dvy.setOptStat("1111");
		fn_rec_kp2_dvz = new F1D("fn_rec_kp2_dvz", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_kp2_dvz.setLineWidth(2);
		fn_rec_kp2_dvz.setLineColor(2);
		fn_rec_kp2_dvz.setOptStat("1111");


		hi_rec_kp2_dp_p =new H2F("hi_rec_kp2_dp_p","hi_rec_kp2_dp_p", 100, 0, 8, 100, -0.05, 0.05);
		hi_rec_kp2_dp_p.setTitleX("p (GeV)");
		hi_rec_kp2_dp_p.setTitleY("#Delta p/p");
		hi_rec_kp2_dp_theta =new H2F("hi_rec_kp2_dp_theta","hi_rec_kp2_dp_theta", 100, 0, 100, 100, -0.05, 0.05);
		hi_rec_kp2_dp_theta.setTitleX("#theta (^o)");
		hi_rec_kp2_dp_theta.setTitleY("#Delta p/p");
		hi_rec_kp2_dp_phi =new H2F("hi_rec_kp2_dp_phi","hi_rec_kp2_dp_phi", 100, 180, -180, 100, -0.05, 0.05);
		hi_rec_kp2_dp_phi.setTitleX("#phi (^o)");
		hi_rec_kp2_dp_phi.setTitleY("#Delta p/p");
		hi_rec_kp2_dp_vz =new H2F("hi_rec_kp2_dp_vz","hi_rec_kp2_dp_vz", 100, -10, 10, 100, -0.05, 0.05);
		hi_rec_kp2_dp_vz.setTitleX("vz");
		hi_rec_kp2_dp_vz.setTitleY("#Delta p/p");
		hi_rec_kp2_dtheta_p =new H2F("hi_rec_kp2_dtheta_p","hi_rec_kp2_dtheta_p", 100, 0, 8, 100, -0.5, 0.5);
		hi_rec_kp2_dtheta_p.setTitleX("p (GeV)");
		hi_rec_kp2_dtheta_p.setTitleY("#Delta#theta (^o)");
		hi_rec_kp2_dtheta_theta =new H2F("hi_rec_kp2_dtheta_theta","hi_rec_kp2_dtheta_theta", 100, 0, 100, 100, -0.5, 0.5);
		hi_rec_kp2_dtheta_theta.setTitleX("#theta (^o)");
		hi_rec_kp2_dtheta_theta.setTitleY("#Delta#theta (^o)");
		hi_rec_kp2_dtheta_phi =new H2F("hi_rec_kp2_dtheta_phi","hi_rec_kp2_dtheta_phi", 100, 180, -180, 100, -0.5, 0.5);
		hi_rec_kp2_dtheta_phi.setTitleX("#phi (^o)");
		hi_rec_kp2_dtheta_phi.setTitleY("#Delta#theta (^o)");
		hi_rec_kp2_dtheta_vz =new H2F("hi_rec_kp2_dtheta_vz","hi_rec_kp2_dtheta_vz", 100, -10, 10, 100, -0.5, 0.5);
		hi_rec_kp2_dtheta_vz.setTitleX("vz");
		hi_rec_kp2_dtheta_vz.setTitleY("#Delta#theta (^o)");
		hi_rec_kp2_dphi_p =new H2F("hi_rec_kp2_dphi_p","hi_rec_kp2_dphi_p", 100, 0, 8, 100, -2.0, 2.0);
		hi_rec_kp2_dphi_p.setTitleX("p (GeV)");
		hi_rec_kp2_dphi_p.setTitleY("#Delta#phi (^o)");
		hi_rec_kp2_dphi_theta =new H2F("hi_rec_kp2_dphi_theta","hi_rec_kp2_dphi_theta", 100, 0, 100, 100, -2.0, 2.0);
		hi_rec_kp2_dphi_theta.setTitleX("#theta (^o)");
		hi_rec_kp2_dphi_theta.setTitleY("#Delta#phi (^o)");
		hi_rec_kp2_dphi_phi =new H2F("hi_rec_kp2_dphi_phi","hi_rec_kp2_dphi_phi", 100, 180, -180, 100, -2.0, 2.0);
		hi_rec_kp2_dphi_phi.setTitleX("#phi (^o)");
		hi_rec_kp2_dphi_phi.setTitleY("#Delta#phi (^o)");
		hi_rec_kp2_dphi_vz =new H2F("hi_rec_kp2_dphi_vz","hi_rec_kp2_dphi_vz", 100, -10, 10, 100, -2.0, 2.0);
		hi_rec_kp2_dphi_vz.setTitleX("vz");
		hi_rec_kp2_dphi_vz.setTitleY("#Delta#phi (^o)");
		hi_rec_kp2_dvz_p =new H2F("hi_rec_kp2_dvz_p","hi_rec_kp2_dvz_p", 100, 0, 8, 100, -3.0, 3.0);
		hi_rec_kp2_dvz_p.setTitleX("p (GeV)");
		hi_rec_kp2_dvz_p.setTitleY("#Delta vz");
		hi_rec_kp2_dvz_theta =new H2F("hi_rec_kp2_dvz_theta","hi_rec_kp2_dvz_theta", 100, 0, 100, 100, -3.0, 3.0);
		hi_rec_kp2_dvz_theta.setTitleX("#theta (^o)");
		hi_rec_kp2_dvz_theta.setTitleY("#Delta vz");
		hi_rec_kp2_dvz_phi =new H2F("hi_rec_kp2_dvz_phi","hi_rec_kp2_dvz_phi", 100, 180, -180, 100, -3.0, 3.0);
		hi_rec_kp2_dvz_phi.setTitleX("#phi (^o)");
		hi_rec_kp2_dvz_phi.setTitleY("#Delta vz");
		hi_rec_kp2_dvz_vz =new H2F("hi_rec_kp2_dvz_vz","hi_rec_kp2_dvz_vz", 100, -10, 10, 100, -3.0, 3.0);
		hi_rec_kp2_dvz_vz.setTitleX("vz");
		hi_rec_kp2_dvz_vz.setTitleY("#Delta vz");



	//KM
		// generated
		hi_mc_km_p_the = new H2F("hi_mc_km_p_the", "hi_mc_km_p_the", 100, 0, 100, 100, 0, 8);
		hi_mc_km_p_the.setTitleX("#theta (^o)");
		hi_mc_km_p_the.setTitleY("p (GeV)");
		hi_mc_km_p_phi = new H2F("hi_mc_km_p_phi", "hi_mc_km_p_phi", 100, -180, 180, 100, 0, 8);
		hi_mc_km_p_phi.setTitleX("#phi (^o)");
		hi_mc_km_p_phi.setTitleY("p (GeV)");
		hi_mc_km_the_phi = new H2F("hi_mc_km_the_phi", "hi_mc_km_the_phi", 100, -180, 180, 100, 0, 100);
		hi_mc_km_the_phi.setTitleX("#phi (^o)");
		hi_mc_km_the_phi.setTitleY("#theta (^o)");

		//reconstructed
		hi_rec_km_p_the = new H2F("hi_rec_km_p_the", "hi_rec_km_p_the", 100, 0, 100, 100, 0, 8);
		hi_rec_km_p_the.setTitleX("#theta (^o)");
		hi_rec_km_p_the.setTitleY("p (GeV)");
		hi_rec_km_p_phi = new H2F("hi_rec_km_p_phi", "hi_rec_km_p_phi", 100, -180, 180, 100, 0, 8);
		hi_rec_km_p_phi.setTitleX("#phi (^o)");
		hi_rec_km_p_phi.setTitleY("p (GeV)");
		hi_rec_km_the_phi = new H2F("hi_rec_km_the_phi", "hi_rec_km_the_phi", 100, -180, 180, 100, 0, 100);
		hi_rec_km_the_phi.setTitleX("#phi (^o)");
		hi_rec_km_the_phi.setTitleY("#theta (^o)");

		hi_rec_km_dvx = new H1F("hi_rec_km_dvx", "hi_rec_km_dvx", 100, -3.0, 3.0);
		hi_rec_km_dvx.setFillColor(LIGHTGREEN);
		hi_rec_km_dvx.setTitleX("#Delta vx");
		hi_rec_km_dvy = new H1F("hi_rec_km_dvy", "hi_rec_km_dvy", 100, -3.0, 3.0);
		hi_rec_km_dvy.setFillColor(LIGHTGREEN);
		hi_rec_km_dvy.setTitleX("#Delta vy");
		hi_rec_km_dvz = new H1F("hi_rec_km_dvz", "hi_rec_km_dvz", 100, -3.0, 3.0);
		hi_rec_km_dvz.setFillColor(LIGHTGREEN);
		hi_rec_km_dvz.setTitleX("#Delta vz");
		hi_rec_km_dp = new H1F("hi_rec_km_dp", "hi_rec_km_dp", 100, -0.05, 0.05);
		hi_rec_km_dp.setTitleX("#delta P/P");
		hi_rec_km_dp.setFillColor(LIGHTGREEN);
		hi_rec_km_dtheta = new H1F("hi_rec_km_dtheta", "hi_rec_km_dtheta", 100, -0.5, 0.5);
		hi_rec_km_dtheta.setTitleX("#Delta#theta (^o)");
		hi_rec_km_dtheta.setFillColor(LIGHTGREEN);
		hi_rec_km_dphi = new H1F("hi_rec_km_dphi","hi_rec_km_dphi", 100, -2.0, 2.0);
		hi_rec_km_dphi.setTitleX("#Delta#phi (^o)");
		hi_rec_km_dphi.setTitleY("count");
		hi_rec_km_dphi.setFillColor(LIGHTGREEN);

		// delta functions 
		fn_rec_km_dp = new F1D("fn_rec_km_dp", "[amp]*gaus(x,[mean],[sigma])", -1.0, 1.0);
		fn_rec_km_dp.setLineWidth(2);
		fn_rec_km_dp.setLineColor(2);
		fn_rec_km_dp.setOptStat("1111");
		fn_rec_km_dtheta = new F1D("fn_rec_km_dtheta", "[amp]*gaus(x,[mean],[sigma])", -0.5, 0.5);
		fn_rec_km_dtheta.setLineWidth(2);
		fn_rec_km_dtheta.setLineColor(2);
		fn_rec_km_dtheta.setOptStat("1111");
		fn_rec_km_dphi = new F1D("fn_rec_km_dphi", "[amp]*gaus(x,[mean],[sigma])", -4, 4);
		fn_rec_km_dphi.setLineWidth(2);
		fn_rec_km_dphi.setLineColor(2);
		fn_rec_km_dphi.setOptStat("1111");
		fn_rec_km_dvx = new F1D("fn_rec_km_dvx", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_km_dvx.setLineWidth(2);
		fn_rec_km_dvx.setLineColor(2);
		fn_rec_km_dvx.setOptStat("1111");
		fn_rec_km_dvy = new F1D("fn_rec_km_dvy", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_km_dvy.setLineWidth(2);
		fn_rec_km_dvy.setLineColor(2);
		fn_rec_km_dvy.setOptStat("1111");
		fn_rec_km_dvz = new F1D("fn_rec_km_dvz", "[amp]*gaus(x,[mean],[sigma])", -5, 5);
		fn_rec_km_dvz.setLineWidth(2);
		fn_rec_km_dvz.setLineColor(2);
		fn_rec_km_dvz.setOptStat("1111");

		hi_rec_km_dp_p =new H2F("hi_rec_km_dp_p","hi_rec_km_dp_p", 100, 0, 8, 100, -0.05, 0.05);
		hi_rec_km_dp_p.setTitleX("p (GeV)");
		hi_rec_km_dp_p.setTitleY("#Delta p/p");
		hi_rec_km_dp_theta =new H2F("hi_rec_km_dp_theta","hi_rec_km_dp_theta", 100, 0, 100, 100, -0.05, 0.05);
		hi_rec_km_dp_theta.setTitleX("#theta (^o)");
		hi_rec_km_dp_theta.setTitleY("#Delta p/p");
		hi_rec_km_dp_phi =new H2F("hi_rec_km_dp_phi","hi_rec_km_dp_phi", 100, 180, -180, 100, -0.05, 0.05);
		hi_rec_km_dp_phi.setTitleX("#phi (^o)");
		hi_rec_km_dp_phi.setTitleY("#Delta p/p");
		hi_rec_km_dp_vz =new H2F("hi_rec_km_dp_vz","hi_rec_km_dp_vz", 100, -10, 10, 100, -0.05, 0.05);
		hi_rec_km_dp_vz.setTitleX("vz");
		hi_rec_km_dp_vz.setTitleY("#Delta p/p");
		hi_rec_km_dtheta_p =new H2F("hi_rec_km_dtheta_p","hi_rec_km_dtheta_p", 100, 0, 8, 100, -0.5, 0.5);
		hi_rec_km_dtheta_p.setTitleX("p (GeV)");
		hi_rec_km_dtheta_p.setTitleY("#Delta#theta (^o)");
		hi_rec_km_dtheta_theta =new H2F("hi_rec_km_dtheta_theta","hi_rec_km_dtheta_theta", 100, 0, 100, 100, -0.5, 0.5);
		hi_rec_km_dtheta_theta.setTitleX("#theta (^o)");
		hi_rec_km_dtheta_theta.setTitleY("#Delta#theta (^o)");
		hi_rec_km_dtheta_phi =new H2F("hi_rec_km_dtheta_phi","hi_rec_km_dtheta_phi", 100, 180, -180, 100, -0.5, 0.5);
		hi_rec_km_dtheta_phi.setTitleX("#phi (^o)");
		hi_rec_km_dtheta_phi.setTitleY("#Delta#theta (^o)");
		hi_rec_km_dtheta_vz =new H2F("hi_rec_km_dtheta_vz","hi_rec_km_dtheta_vz", 100, -10, 10, 100, -0.5, 0.5);
		hi_rec_km_dtheta_vz.setTitleX("vz");
		hi_rec_km_dtheta_vz.setTitleY("#Delta#theta (^o)");
		hi_rec_km_dphi_p =new H2F("hi_rec_km_dphi_p","hi_rec_km_dphi_p", 100, 0, 8, 100, -2.0, 2.0);
		hi_rec_km_dphi_p.setTitleX("p (GeV)");
		hi_rec_km_dphi_p.setTitleY("#Delta#phi (^o)");
		hi_rec_km_dphi_theta =new H2F("hi_rec_km_dphi_theta","hi_rec_km_dphi_theta", 100, 0, 100, 100, -2.0, 2.0);
		hi_rec_km_dphi_theta.setTitleX("#theta (^o)");
		hi_rec_km_dphi_theta.setTitleY("#Delta#phi (^o)");
		hi_rec_km_dphi_phi =new H2F("hi_rec_km_dphi_phi","hi_rec_km_dphi_phi", 100, 180, -180, 100, -2.0, 2.0);
		hi_rec_km_dphi_phi.setTitleX("#phi (^o)");
		hi_rec_km_dphi_phi.setTitleY("#Delta#phi (^o)");
		hi_rec_km_dphi_vz =new H2F("hi_rec_km_dphi_vz","hi_rec_km_dphi_vz", 100, -10, 10, 100, -2.0, 2.0);
		hi_rec_km_dphi_vz.setTitleX("vz");
		hi_rec_km_dphi_vz.setTitleY("#Delta#phi (^o)");
		hi_rec_km_dvz_p =new H2F("hi_rec_km_dvz_p","hi_rec_km_dvz_p", 100, 0, 8, 100, -3.0, 3.0);
		hi_rec_km_dvz_p.setTitleX("p (GeV)");
		hi_rec_km_dvz_p.setTitleY("#Delta vz");
		hi_rec_km_dvz_theta =new H2F("hi_rec_km_dvz_theta","hi_rec_km_dvz_theta", 100, 0, 100, 100, -3.0, 3.0);
		hi_rec_km_dvz_theta.setTitleX("#theta (^o)");
		hi_rec_km_dvz_theta.setTitleY("#Delta vz");
		hi_rec_km_dvz_phi =new H2F("hi_rec_km_dvz_phi","hi_rec_km_dvz_phi", 100, 180, -180, 100, -3.0, 3.0);
		hi_rec_km_dvz_phi.setTitleX("#phi (^o)");
		hi_rec_km_dvz_phi.setTitleY("#Delta vz");
		hi_rec_km_dvz_vz =new H2F("hi_rec_km_dvz_vz","hi_rec_km_dvz_vz", 100, -10, 10, 100, -3.0, 3.0);
		hi_rec_km_dvz_vz.setTitleX("vz");
		hi_rec_km_dvz_vz.setTitleY("#Delta vz");


	//counter for no of rec particles 
		hi_pip_counter = new H1F("hi_pip_counter", "hi_pip_counter", 5, 0, 5);
		hi_pip_counter.setFillColor(LIGHTGREEN);
		hi_pim_counter = new H1F("hi_pim_counter", "hi_pim_counter", 5, 0, 5);
		hi_pim_counter.setFillColor(LIGHTGREEN);
		hi_kp_counter = new H1F("hi_kp_counter", "hi_kp_counter", 5, 0, 5);
		hi_kp_counter.setFillColor(LIGHTGREEN);
		hi_km_counter = new H1F("hi_km_counter", "hi_km_counter", 5, 0, 5);
		hi_km_counter.setFillColor(LIGHTGREEN);
		hi_prot_counter = new H1F("hi_prot_counter", "hi_prot_counter", 5, 0, 5);
		hi_prot_counter.setFillColor(LIGHTGREEN);
		hi_fpip_counter = new H1F("hi_fpip_counter", "hi_fpip_counter", 5, 0, 5);
		hi_fpip_counter.setFillColor(LIGHTGREEN);
		hi_fpim_counter = new H1F("hi_fpim_counter", "hi_fpim_counter", 5, 0, 5);
		hi_fpim_counter.setFillColor(LIGHTGREEN);
		hi_fkp_counter = new H1F("hi_fkp_counter", "hi_fkp_counter", 5, 0, 5);
		hi_fkp_counter.setFillColor(LIGHTGREEN);
		hi_fkm_counter = new H1F("hi_fkm_counter", "hi_fkm_counter", 5, 0, 5);
		hi_fkm_counter.setFillColor(LIGHTGREEN);
		hi_fprot_counter = new H1F("hi_fprot_counter", "hi_fprot_counter", 5, 0, 5);
		hi_fprot_counter.setFillColor(LIGHTGREEN);
		hi_cpip_counter = new H1F("hi_cpip_counter", "hi_cpip_counter", 5, 0, 5);
		hi_cpip_counter.setFillColor(LIGHTGREEN);
		hi_cpim_counter = new H1F("hi_cpim_counter", "hi_cpim_counter", 5, 0, 5);
		hi_cpim_counter.setFillColor(LIGHTGREEN);
		hi_ckp_counter = new H1F("hi_ckp_counter", "hi_ckp_counter", 5, 0, 5);
		hi_ckp_counter.setFillColor(LIGHTGREEN);
		hi_ckm_counter = new H1F("hi_ckm_counter", "hi_ckm_counter", 5, 0, 5);
		hi_ckm_counter.setFillColor(LIGHTGREEN);
		hi_cprot_counter = new H1F("hi_cprot_counter", "hi_cprot_counter", 5, 0, 5);
		hi_cprot_counter.setFillColor(LIGHTGREEN);
		
	//ekpkpreqkm	reconstructed
		
		hi_ekpkpkm_MM_ekpkp = new H1F("hi_ekpkpkm_MM_ekpkp", "hi_ekpkpkm_MM_ekpkp", 50, 1.6, 2.1);
		hi_ekpkpkm_MM_ekpkp.setTitle("MM");
		hi_ekpkpkm_MM_ekpkp.setTitleX("MM(eK^+K^+) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkp.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkp.setFillColor(LIGHTGREEN);

		f1_xi = new F1D("f1_xi", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		f1_xi.setParameter(0, 0);
    	f1_xi.setParameter(1, 1);
    	f1_xi.setParameter(2, 0.2);
    	f1_xi.setLineWidth(2);
    	f1_xi.setLineColor(2);
    	f1_xi.setOptStat("1111");

    	// for one with common vertex correction for electron
    	hi_ekpkpkm_MM_ekpkp_withcorrn = new H1F("hi_ekpkpkm_MM_ekpkp_withcorrn", "hi_ekpkpkm_MM_ekpkp_withcorrn", 50, 1.6, 2.1);
		hi_ekpkpkm_MM_ekpkp_withcorrn.setTitle("MM");
		hi_ekpkpkm_MM_ekpkp_withcorrn.setTitleX("MM(eK^+K^+) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkp_withcorrn.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkp_withcorrn.setFillColor(LIGHTGREEN);

		fn_xi_with_cor = new F1D("fn_xi_with_cor", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		fn_xi_with_cor.setParameter(0, 0);
    	fn_xi_with_cor.setParameter(1, 1);
    	fn_xi_with_cor.setParameter(2, 0.2);
    	fn_xi_with_cor.setLineWidth(2);
    	fn_xi_with_cor.setLineColor(2);
    	fn_xi_with_cor.setOptStat("1111");

    //ekpkpkm generated
    	
    	hi_mc_ekpkpkm_mm_ekpkp = new H1F("hi_mc_ekpkpkm_mm_ekpkp", "hi_mc_ekpkpkm_mm_ekpkp", 50, 1.6, 2.1);
		hi_mc_ekpkpkm_mm_ekpkp.setTitle("MM");
		hi_mc_ekpkpkm_mm_ekpkp.setTitleX("MM(eK^+K^+) [GeV/c^2]");
		hi_mc_ekpkpkm_mm_ekpkp.setTitleY("Events/[10 MeV/c^2]");
		hi_mc_ekpkpkm_mm_ekpkp.setFillColor(LIGHTGREEN);

		f1_mc_xi = new F1D("f1_mc_xi", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		//f1_xi.setParameter(0, 0);
    	//f1_xi.setParameter(1, 1);
    	//f1_xi.setParameter(2, 0.2);
    	f1_mc_xi.setLineWidth(2);
    	f1_mc_xi.setLineColor(2);
    	f1_mc_xi.setOptStat("1111");

    	hi_mc_ekpkpkm_mm_ekpkpkm = new H1F("hi_mc_ekpkpkm_mm_ekpkpkm", "hi_mc_ekpkpkm_mm_ekpkpkm", 50, 0.9, 1.4);
		hi_mc_ekpkpkm_mm_ekpkpkm.setTitle("MM");
		hi_mc_ekpkpkm_mm_ekpkpkm.setTitleX("MM(eK^+K^+K^-) [GeV/c^2]");
		hi_mc_ekpkpkm_mm_ekpkpkm.setTitleY("Events/[8 MeV/c^2]");
		hi_mc_ekpkpkm_mm_ekpkpkm.setFillColor(LIGHTGREEN);

    	L_lamda = new DataLine(1.115683, 0, 1.115683, 70);
    	L_lamda.setLineColor(2);
        L_lamda.setLineWidth(2);
        L_lamda.setArrowSizeOrigin(0);
        L_lamda.setArrowSizeEnd(0);
        L_lamda.setArrowAngle(25);
	//	L_lamda.setOrigin(1.115683, 0);
	//	L_lamda.setEnd(1.115683, 700);
	//	L_lamda.setLineColor(2);
		
		L_sigma = new DataLine(1.192642, 0, 1.192642, 70);
		L_sigma.setLineColor(2);
        L_sigma.setLineWidth(2);
        L_sigma.setArrowSizeOrigin(0);
        L_sigma.setArrowSizeEnd(0);
        L_sigma.setArrowAngle(25);
	//	L_sigma.setOrigin(1.192642, 0); //1.18937
	//	L_sigma.setEnd(1.192642, 700); //1.18937
	//	L_sigma.setLineColor(3);

// ekpkpkm MM reconstruction
		hi_ekpkpkm_MM_ekpkpkm = new H1F("hi_ekpkpkm_MM_ekpkpkm", "hi_ekpkpkm_MM_ekpkpkm", 50, 0.9, 1.4);
		hi_ekpkpkm_MM_ekpkpkm.setTitle("MM");
		hi_ekpkpkm_MM_ekpkpkm.setTitleX("MM(eK^+K^+K^-) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm.setTitleY("Events/[8 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm.setFillColor(LIGHTGREEN);

		// for one with common vertex correction for electron
		hi_ekpkpkm_MM_ekpkpkm_withcorrn = new H1F("hi_ekpkpkm_MM_ekpkpkm_withcorrn", "hi_ekpkpkm_MM_ekpkpkm_withcorrn", 50, 0.9, 1.4);
		hi_ekpkpkm_MM_ekpkpkm_withcorrn.setTitle("MM");
		hi_ekpkpkm_MM_ekpkpkm_withcorrn.setTitleX("MM(eK^+K^+K^-) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm_withcorrn.setTitleY("Events/[8 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm_withcorrn.setFillColor(LIGHTGREEN);

		hi_ekpkpkm_IM_kmlambda = new H1F("hi_ekpkpkm_IM_kmlambda", "hi_ekpkpkm_IM_kmlambda", 50, 1.6, 2.1);
		hi_ekpkpkm_IM_kmlambda.setTitle("M");
		hi_ekpkpkm_IM_kmlambda.setTitleX("M(K^-Lambda) [GeV/c^2]");
		hi_ekpkpkm_IM_kmlambda.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_IM_kmlambda.setFillColor(LIGHTGREEN);

		fn_im_kmlambda = new F1D("fn_im_kmlambda", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		//f1_xi.setParameter(0, 0);
    	//f1_xi.setParameter(1, 1);
    	//f1_xi.setParameter(2, 0.2);
    	fn_im_kmlambda.setLineWidth(2);
    	fn_im_kmlambda.setLineColor(2);
    	fn_im_kmlambda.setOptStat("1111");

		hi_ekpkpkm_IM_kmsigma = new H1F("hi_ekpkpkm_IM_kmsigma", "hi_ekpkpkm_IM_kmsigma", 50, 1.6, 2.1);
		hi_ekpkpkm_IM_kmsigma.setTitle("M");
		hi_ekpkpkm_IM_kmsigma.setTitleX("M(K^-Sigma) [GeV/c^2]");
		hi_ekpkpkm_IM_kmsigma.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_IM_kmsigma.setFillColor(LIGHTGREEN);

		fn_im_kmsigma = new F1D("fn_im_kmsigma", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		//f1_xi.setParameter(0, 0);
    	//f1_xi.setParameter(1, 1);
    	//f1_xi.setParameter(2, 0.2);
    	fn_im_kmsigma.setLineWidth(2);
    	fn_im_kmsigma.setLineColor(2);
    	fn_im_kmsigma.setOptStat("1111");


		f1_lambda = new F1D("f1_lambda", "[amp]*gaus(x,[mean],[sigma])", 1.04, 1.14);
		f1_lambda.setParameter(0, 0);
    	f1_lambda.setParameter(1, 1);
    	f1_lambda.setParameter(2, 0.2);
    	f1_lambda.setLineWidth(2);
    	f1_lambda.setLineColor(2);
    	f1_lambda.setOptStat("1111");

    	f1_sigma = new F1D("f1_sigma", "[amp]*gaus(x,[mean],[sigma])", 1.17, 1.24);
		f1_sigma.setParameter(0, 0);
    	f1_sigma.setParameter(1, 1);
    	f1_sigma.setParameter(2, 0.2);
    	f1_sigma.setLineWidth(2);
    	f1_sigma.setLineColor(2);
    	f1_sigma.setOptStat("1111");

    	f_gaus = new F1D("f_gaus", "[amp]*gaus(x,[mean],[sigma])", -2.5, 2.5);
    	f_gaus.setParameter(0, 0);
    	f_gaus.setParameter(1, 1);
    	f_gaus.setParameter(2, 1);
    	f_gaus.setLineWidth(2);
    	f_gaus.setLineColor(2);
    	f_gaus.setOptStat("1111");

    	// ekpkpkmprot detected
    	hi_ekpkpkmprot_MM2 = new H1F("hi_ekpkpkmprot_MM2", "hi_ekpkpkmprot_MM2", 200, -0.5, 0.5);
		hi_ekpkpkmprot_MM2.setTitle("MM2");
		hi_ekpkpkmprot_MM2.setTitleX("MM2(eK^+K^+K^-p) [GeV^2/c^4]");
		hi_ekpkpkmprot_MM2.setTitleY("Events");
		hi_ekpkpkmprot_MM2.setFillColor(LIGHTGREEN);
		
		H_efkpckp_MM_req_km = new H1F("H_efkpckp_MM_req_km", "H_efkpckp_MM_req_km", 50, 1.6, 2.1);
		H_efkpckp_MM_req_km.setTitle("MM requiring an additional K^- anywhere");
		H_efkpckp_MM_req_km.setTitleX("MM(efk^+ck^+) (GeV)");
		H_efkpckp_MM_req_km.setTitleY("count");
		H_efkpckp_MM_req_km.setFillColor(LIGHTGREEN);
		
		H_efkpckp_MM_req_ckm = new H1F("H_efkpckp_MM_req_ckm", "H_efkpckp_MM_req_ckm", 50, 1.6, 2.1);
		H_efkpckp_MM_req_ckm.setTitle("MM requiring an additional K^- in CD");
		H_efkpckp_MM_req_ckm.setTitleX("MM(efk^+ck^+) (GeV)");
		H_efkpckp_MM_req_ckm.setTitleY("count");
		H_efkpckp_MM_req_ckm.setFillColor(LIGHTGREEN);

		H_ekpkp_MM_req_pim = new H1F("H_ekpkp_MM_req_pim", "H_ekpkp_MM_req_pim", 50, 1.6, 2.1);
		H_ekpkp_MM_req_pim.setTitle("MM requiring an additional #pi^-");
		H_ekpkp_MM_req_pim.setTitleX("MM(ek^+k^+) (GeV)");
		H_ekpkp_MM_req_pim.setTitleY("count");
		H_ekpkp_MM_req_pim.setFillColor(LIGHTGREEN);
		
		H_efkpckp_MM_req_pim = new H1F("H_efkpckp_MM_req_pim", "H_efkpckp_MM_req_pim", 50, 1.6, 2.1);
		H_efkpckp_MM_req_pim.setTitle("MM requiring an additional #pi^- anywhere");
		H_efkpckp_MM_req_pim.setTitleX("MM(efk^+ck^+) (GeV)");
		H_efkpckp_MM_req_pim.setTitleY("count");
		H_efkpckp_MM_req_pim.setFillColor(LIGHTGREEN);
		
		H_efkpckp_MM_req_cpim = new H1F("H_efkpckp_MM_req_cpim", "H_efkpckp_MM_req_cpim", 50, 1.6, 2.1);
		H_efkpckp_MM_req_cpim.setTitle("MM requiring an additional #pi^- in CD");
		H_efkpckp_MM_req_cpim.setTitleX("MM(efk^+ck^+) (GeV)");
		H_efkpckp_MM_req_cpim.setTitleY("count");
		H_efkpckp_MM_req_cpim.setFillColor(LIGHTGREEN);
		
		//FD particle vertex
		H_pim_vz = new H1F("H_pim_vz", "H_pim_vz", 100, -20, 20);
		H_pim_vz.setFillColor(LIGHTGREEN);
		H_pip_vz = new H1F("H_pip_vz", "H_pip_vz", 100, -20, 20);
		H_pip_vz.setFillColor(LIGHTGREEN);
		H_kp_vz = new H1F("H_kp_vz", "H_kp_vz", 100, -20, 20);
		H_kp_vz.setFillColor(LIGHTGREEN);
		H_km_vz = new H1F("H_km_vz", "H_km_vz", 100, -20, 20);
		H_km_vz.setFillColor(LIGHTGREEN);
		H_prot_vz = new H1F("H_prot_vz", "H_prot_vz", 100, -20, 20);
		H_prot_vz.setFillColor(LIGHTGREEN);
		
		//FD particle path
		H_prot_FTOF1b_path = new H1F("H_prot_FTOF1b_path", "H_prot_FTOF1b_path", 200, 400, 900);
		H_prot_FTOF1b_path.setFillColor(LIGHTGREEN);
		H_pip_FTOF1b_path = new H1F("H_pip_FTOF1b_path", "H_pip_FTOF1b_path", 200, 400, 900);
		H_pip_FTOF1b_path.setFillColor(LIGHTGREEN);
		H_pim_FTOF1b_path = new H1F("H_pim_FTOF1b_path", "H_pim_FTOF1b_path", 200, 400, 900);
		H_pim_FTOF1b_path.setFillColor(LIGHTGREEN);
		H_kp_FTOF1b_path = new H1F("H_kp_FTOF1b_path", "H_kp_FTOF1b_path", 200, 400, 900);
		H_kp_FTOF1b_path.setFillColor(LIGHTGREEN);
		H_km_FTOF1b_path = new H1F("H_km_FTOF1b_path", "H_km_FTOF1b_path", 200, 400, 900);
		H_km_FTOF1b_path.setFillColor(LIGHTGREEN);
		//FD particle time
		H_prot_FTOF1b_t = new H1F("H_prot_FTOF1b_t", "H_prot_FTOF1b_t", 200, 0, 80);
		H_prot_FTOF1b_t.setFillColor(LIGHTGREEN);
		H_pip_FTOF1b_t = new H1F("H_pip_FTOF1b_t", "H_pip_FTOF1b_t", 200, 0, 80);
		H_pip_FTOF1b_t.setFillColor(LIGHTGREEN);
		H_pim_FTOF1b_t = new H1F("H_pim_FTOF1b_t", "H_pim_FTOF1b_t", 200, 0, 80);
		H_pim_FTOF1b_t.setFillColor(LIGHTGREEN);
		H_kp_FTOF1b_t = new H1F("H_kp_FTOF1b_t", "H_kp_FTOF1b_t", 200, 0, 80);
		H_kp_FTOF1b_t.setFillColor(LIGHTGREEN);
		H_km_FTOF1b_t = new H1F("H_km_FTOF1b_t", "H_km_FTOF1b_t", 200, 0, 80);
		H_km_FTOF1b_t.setFillColor(LIGHTGREEN);
		
		// FTOF particles vetrex time vs momentum
		H_pip_vt_p = new H2F("H_pip_vt_p","H_pip_vt_p", 100, -0.5, 0.5, 100, 0, 10.6);
		H_pip_vt_p.setTitle("pip vt vs mom");
		H_pip_vt_p.setTitleX("vt (ns)");
		H_pip_vt_p.setTitleY("p (GeV)");
		H_pim_vt_p = new H2F("H_pim_vt_p","H_pim_vt_p", 100, -0.5, 0.5, 100, 0, 10.6);
		H_pim_vt_p.setTitle("pim vt vs mom");
		H_pim_vt_p.setTitleX("vt (ns)");
		H_pim_vt_p.setTitleY("p (GeV)");
		H_kp_vt_p = new H2F("H_kp_vt_p","H_kp_vt_p", 100, -0.5, 0.5, 100, 0, 10.6);
		H_kp_vt_p.setTitle("kp vt vs mom");
		H_kp_vt_p.setTitleX("vt (ns)");
		H_kp_vt_p.setTitleY("p (GeV)");
		H_km_vt_p = new H2F("H_km_vt_p","H_km_vt_p", 100, -0.5, 0.5, 100, 0, 10.6);
		H_km_vt_p.setTitle("km vt vs mom");
		H_km_vt_p.setTitleX("vt (ns)");
		H_km_vt_p.setTitleY("p (GeV)");
		H_prot_vt_p = new H2F("H_prot_vt_p","H_prot_vt_p", 100, -0.5, 0.5, 100, 0, 10.6);
		H_prot_vt_p.setTitle("prot vt vs mom");
		H_prot_vt_p.setTitleX("vt (ns)");
		H_prot_vt_p.setTitleY("p (GeV)");
		
		//CD particle vertex
		H_pimc_vz = new H1F("H_pimc_vz", "H_pimc_vz", 100, -20, 20);
		H_pimc_vz.setFillColor(LIGHTGREEN);
		H_pipc_vz = new H1F("H_pipc_vz", "H_pipc_vz", 100, -20, 20);
		H_pipc_vz.setFillColor(LIGHTGREEN);
		H_kpc_vz = new H1F("H_kpc_vz", "H_kpc_vz", 100, -20, 20);
		H_kpc_vz.setFillColor(LIGHTGREEN);
		H_kmc_vz = new H1F("H_km_vz", "H_kmc_vz", 100, -20, 20);
		H_kmc_vz.setFillColor(LIGHTGREEN);
		H_protc_vz = new H1F("H_protc_vz", "H_protc_vz", 100, -20, 20);
		H_protc_vz.setFillColor(LIGHTGREEN);
		
		//CD particle path
		H_prot_CTOF_path = new H1F("H_prot_CTOF_path", "H_prot_CTOF_path", 200, 0, 100);
		H_prot_CTOF_path.setFillColor(LIGHTGREEN);
		H_pip_CTOF_path = new H1F("H_pip_CTOF_path", "H_pip_CTOF_path", 200, 0, 100);
		H_pip_CTOF_path.setFillColor(LIGHTGREEN);
		H_pim_CTOF_path = new H1F("H_pim_CTOF_path", "H_pim_CTOF_path", 200, 0, 100);
		H_pim_CTOF_path.setFillColor(LIGHTGREEN);
		H_kp_CTOF_path = new H1F("H_kp_CTOF_path", "H_kp_CTOF_path", 200, 0, 100);
		H_kp_CTOF_path.setFillColor(LIGHTGREEN);
		H_km_CTOF_path = new H1F("H_km_CTOF_path", "H_km_CTOF_path", 200, 0, 100);
		H_km_CTOF_path.setFillColor(LIGHTGREEN);
		
		//CD particle time
		H_prot_CTOF_t = new H1F("H_prot_CTOF_t", "H_prot_CTOF_t", 200, -5, 20);
		H_prot_CTOF_t.setFillColor(LIGHTGREEN);
		H_pip_CTOF_t = new H1F("H_pip_CTOF_t", "H_pip_CTOF_t", 200, -5, 20);
		H_pip_CTOF_t.setFillColor(LIGHTGREEN);
		H_pim_CTOF_t = new H1F("H_pim_CTOF_t", "H_pim_CTOF_t", 200, -5, 20);
		H_pim_CTOF_t.setFillColor(LIGHTGREEN);
		H_kp_CTOF_t = new H1F("H_kp_CTOF_t", "H_kp_CTOF_t", 200, -5, 20);
		H_kp_CTOF_t.setFillColor(LIGHTGREEN);
		H_km_CTOF_t = new H1F("H_km_CTOF_t", "H_km_CTOF_t", 200, -5, 20);
		H_km_CTOF_t.setFillColor(LIGHTGREEN);		
		
		
		//CTOF particles vetrex time vs momentum
		H_pipc_vt_p = new H2F("H_pipc_vt_p","H_pipc_vt_p", 100, -1, 1, 100, 0, 4);
		H_pipc_vt_p.setTitle("pipc vt vs mom");
		H_pipc_vt_p.setTitleX("vt (ns)");
		H_pipc_vt_p.setTitleY("p (GeV)");
		
		H_pimc_vt_p = new H2F("H_pimc_vt_p","H_pimc_vt_p", 100, -1, 1, 100, 0, 4);
		H_pimc_vt_p.setTitle("pimc vt vs mom");
		H_pimc_vt_p.setTitleX("vt (ns)");
		H_pimc_vt_p.setTitleY("p (GeV)");
		
		H_kpc_vt_p = new H2F("H_kpc_vt_p","H_kpc_vt_p", 100, -1, 1, 100, 0, 4);
		H_kpc_vt_p.setTitle("kpc vt vs mom");
		H_kpc_vt_p.setTitleX("vt (ns)");
		H_kpc_vt_p.setTitleY("p (GeV)");
		
		H_kmc_vt_p = new H2F("H_kmc_vt_p","H_kmc_vt_p", 100, -1, 1, 100, 0, 4);
		H_kmc_vt_p.setTitle("kmc vt vs mom");
		H_kmc_vt_p.setTitleX("vt (ns)");
		H_kmc_vt_p.setTitleY("p (GeV)");
		
		H_protc_vt_p = new H2F("H_protc_vt_p","H_protc_vt_p", 100, -1, 1, 100, 0, 4);
		H_protc_vt_p.setTitle("protc vt vs mom");
		H_protc_vt_p.setTitleX("vt (ns)");
		H_protc_vt_p.setTitleY("p (GeV)");
		
		
		// FD particle beta vs momentum by charge
		H_FD_pos_beta_mom = new H2F("H_FD_pos_beta_mom", "H_FD_pos_beta_mom", 100, 0, 10.6, 100, 0, 1.2);
		H_FD_pos_beta_mom.setTitle("POS  #beta vs mom");
		H_FD_pos_beta_mom.setTitleX("p (GeV)");
		H_FD_pos_beta_mom.setTitleY("FTB #beta");
		H_FD_neg_beta_mom = new H2F("H_FD_neg_beta_mom", "H_FD_neg_beta_mom", 100, 0, 10.6, 100, 0, 1.2);
		H_FD_neg_beta_mom.setTitle("NEG  #beta vs mom");
		H_FD_neg_beta_mom.setTitleX("p (GeV)");
		H_FD_neg_beta_mom.setTitleY("FTB #beta");
		H_FD_neutral_beta_mom = new H2F("H_FD_neutral_beta_mom", "H_FD_neutral_beta_mom", 100, 0, 10.6, 100, 0, 1.2);
		H_FD_neutral_beta_mom.setTitle("NEUTRAL  #beta vs mom");
		H_FD_neutral_beta_mom.setTitleX("p (GeV)");
		H_FD_neutral_beta_mom.setTitleY("FTB #beta");
		hi_fd_pos_mass = new H1F("hi_fd_pos_mass", "hi_fd_pos_mass", 150, -0.5, 4.5);
		hi_fd_pos_mass.setFillColor(LIGHTGREEN);
		H_FD_pos_mass_mom = new H2F("H_FD_pos_mass_mom", "H_FD_pos_mass_mom", 100, 0, 7, 150, -0.5, 4.5);
		H_FD_pos_mass_mom.setTitle("POS Mass^2 vs mom");
		H_FD_pos_mass_mom.setTitleX("p (GeV)");
		H_FD_pos_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_fd_neg_mass = new H1F("hi_fd_neg_mass", "hi_fd_neg_mass", 150, -0.5, 2.0);
		hi_fd_neg_mass.setFillColor(LIGHTGREEN);
		H_FD_neg_mass_mom = new H2F("H_FD_neg_mass_mom", "H_FD_neg_mass_mom", 100, 0, 7, 150, -0.5, 2.0);
		H_FD_neg_mass_mom.setTitle("NEG Mass^2 vs mom");
		H_FD_neg_mass_mom.setTitleX("p (GeV)");
		H_FD_neg_mass_mom.setTitleY("M^2 (GeV^2)");
		H_FD_neutral_mass_mom = new H2F("H_FD_neutral_mass_mom", "H_FD_neutral_mass_mom", 100, 0, 7, 150, -0.5, 2.0);
		H_FD_neutral_mass_mom.setTitle("NEUTRAL Mass^2 vs mom");
		H_FD_neutral_mass_mom.setTitleX("p (GeV)");
		H_FD_neutral_mass_mom.setTitleY("M^2 (GeV^2)");
		H_FD_pos_mass_the = new H2F("H_FD_pos_mass_the", "H_FD_pos_mass_the", 100, 0, 80, 100, -0.5, 4.5);
		H_FD_pos_mass_the.setTitle("POS Mass^2 vs #theta");
		H_FD_pos_mass_the.setTitleX("#theta (^o)");
		H_FD_pos_mass_the.setTitleY("M^2 (GeV^2");
		H_FD_neg_mass_the = new H2F("H_FD_neg_mass_the", "H_FD_neg_mass_the", 100, 0, 80, 100, -0.5, 2);
		H_FD_neg_mass_the.setTitle("NEG Mass^2 vs #theta");
		H_FD_neg_mass_the.setTitleX("#theta (^o)");
		H_FD_neg_mass_the.setTitleY("M^2 (GeV^2");
		H_FD_neutral_mass_the = new H2F("H_FD_neutral_mass_the", "H_FD_neutral_mass_the", 100, 0, 80, 100, -0.5, 2);
		H_FD_neutral_mass_the.setTitle("NEUTRAL Mass^2 vs #theta");
		H_FD_neutral_mass_the.setTitleX("#theta (^o)");
		H_FD_neutral_mass_the.setTitleY("M^2 (GeV^2");
		H_FD_pos_mass_phi = new H2F("H_FD_pos_mass_phi", "H_FD_pos_mass_phi", 100, -180, 180, 100, -0.5, 4.5);
		H_FD_pos_mass_phi.setTitle("POS Mass^2 vs #phi");
		H_FD_pos_mass_phi.setTitleX("#phi (^o)");
		H_FD_pos_mass_phi.setTitleY("M^2 (GeV^2");
		H_FD_neg_mass_phi = new H2F("H_FD_neg_mass_phi", "H_FD_neg_mass_phi", 100, -180, 180, 100, -0.5, 2);
		H_FD_neg_mass_phi.setTitle("NEG Mass^2 vs #phi");
		H_FD_neg_mass_phi.setTitleX("#phi (^o)");
		H_FD_neg_mass_phi.setTitleY("M^2 (GeV^2");
		H_FD_neutral_mass_phi = new H2F("H_FD_neutral_mass_phi", "H_FD_neutral_mass_phi", 100, -180, 180, 100, -0.5, 2);
		H_FD_neutral_mass_phi.setTitle("NEUTRAL Mass^2 vs #phi");
		H_FD_neutral_mass_phi.setTitleX("#phi (^o)");
		H_FD_neutral_mass_phi.setTitleY("M^2 (GeV^2");
		
		
		// CD particle beta vs mom by charge
		H_CD_pos_beta_mom = new H2F("H_CD_pos_beta_mom", "H_CD_pos_beta_mom", 100, 0, 4.0, 100, 0, 1.2);
		H_CD_pos_beta_mom.setTitle("POS  #beta vs mom");
		H_CD_pos_beta_mom.setTitleX("p (GeV)");
		H_CD_pos_beta_mom.setTitleY("FTB #beta");
		H_CD_neg_beta_mom = new H2F("H_CD_neg_beta_mom", "H_CD_neg_beta_mom", 100, 0, 4.0, 100, 0, 1.2);
		H_CD_neg_beta_mom.setTitle("NEG  #beta vs mom");
		H_CD_neg_beta_mom.setTitleX("p (GeV)");
		H_CD_neg_beta_mom.setTitleY("FTB #beta");
		H_CD_neutral_beta_mom = new H2F("H_CD_neutral_beta_mom", "H_CD_neutral_beta_mom", 100, 0, 4.0, 100, 0, 1.2);
		H_CD_neutral_beta_mom.setTitle("NEUTRAL  #beta vs mom");
		H_CD_neutral_beta_mom.setTitleX("p (GeV)");
		H_CD_neutral_beta_mom.setTitleY("FTB #beta");
		hi_cd_pos_mass = new H1F("hi_cd_pos_mass", "hi_cd_pos_mass", 150, -0.5, 4.5);
		hi_cd_pos_mass.setFillColor(LIGHTGREEN);
		H_CD_pos_mass_mom = new H2F("H_CD_pos_mass_mom", "H_CD_pos_mass_mom", 100, 0, 4, 150, -0.5, 4.5);
		H_CD_pos_mass_mom.setTitle("POS Mass^2 vs mom");
		H_CD_pos_mass_mom.setTitleX("p (GeV)");
		H_CD_pos_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_cd_neg_mass = new H1F("hi_fd_neg_mass", "hi_fd_neg_mass", 150, -0.5, 2.0);
		hi_cd_neg_mass.setFillColor(LIGHTGREEN);
		H_CD_neg_mass_mom = new H2F("H_CD_neg_mass_mom", "H_CD_neg_mass_mom", 100, 0, 4, 150, -0.5, 2.0);
		H_CD_neg_mass_mom.setTitle("NEG Mass^2 vs mom");
		H_CD_neg_mass_mom.setTitleX("p (GeV)");
		H_CD_neg_mass_mom.setTitleY("M^2 (GeV^2)");
		H_CD_neutral_mass_mom = new H2F("H_CD_neutral_mass_mom", "H_CD_neutral_mass_mom", 100, 0, 4, 150, -0.5, 2.0);
		H_CD_neutral_mass_mom.setTitle("NEUTRAL Mass^2 vs mom");
		H_CD_neutral_mass_mom.setTitleX("p (GeV)");
		H_CD_neutral_mass_mom.setTitleY("M^2 (GeV^2)");
		H_CD_pos_mass_the = new H2F("H_CD_pos_mass_the", "H_CD_pos_mass_the", 150, 30, 180, 100, -0.5, 4.5);
		H_CD_pos_mass_the.setTitle("POS Mass^2 vs #theta");
		H_CD_pos_mass_the.setTitleX("#theta (^o)");
		H_CD_pos_mass_the.setTitleY("M^2 (GeV^2");
		H_CD_neg_mass_the = new H2F("H_CD_neg_mass_the", "H_CD_neg_mass_the", 150, 30, 180, 100, -0.5, 2);
		H_CD_neg_mass_the.setTitle("NEG Mass^2 vs #theta");
		H_CD_neg_mass_the.setTitleX("#theta (^o)");
		H_CD_neg_mass_the.setTitleY("M^2 (GeV^2");
		H_CD_neutral_mass_the = new H2F("H_CD_neutral_mass_the", "H_CD_neutral_mass_the", 150, 30, 180, 100, -0.5, 2);
		H_CD_neutral_mass_the.setTitle("NEUTRAL Mass^2 vs #theta");
		H_CD_neutral_mass_the.setTitleX("#theta (^o)");
		H_CD_neutral_mass_the.setTitleY("M^2 (GeV^2");
		H_CD_pos_mass_phi = new H2F("H_CD_pos_mass_phi", "H_CD_pos_mass_phi", 150, -180, 180, 100, -0.5, 4.5);
		H_CD_pos_mass_phi.setTitle("POS Mass^2 vs #phi");
		H_CD_pos_mass_phi.setTitleX("#phi (^o)");
		H_CD_pos_mass_phi.setTitleY("M^2 (GeV^2");
		H_CD_neg_mass_phi = new H2F("H_CD_neg_mass_phi", "H_CD_neg_mass_phi", 150, -180, 180, 100, -0.5, 2);
		H_CD_neg_mass_phi.setTitle("NEG Mass^2 vs #phi");
		H_CD_neg_mass_phi.setTitleX("#phi (^o)");
		H_CD_neg_mass_phi.setTitleY("M^2 (GeV^2");
		H_CD_neutral_mass_phi = new H2F("H_CD_neutral_mass_phi", "H_CD_neutral_mass_phi", 150, -180, 180, 100, -0.5, 2);
		H_CD_neutral_mass_phi.setTitle("NEUTRAL Mass^2 vs #phi");
		H_CD_neutral_mass_phi.setTitleX("#phi (^o)");
		H_CD_neutral_mass_phi.setTitleY("M^2 (GeV^2");
		

		// rec electron
		dg_rec_electron = new DataGroup(3, 5);
		dg_rec_electron.addDataSet(H_FT_e_p_the, 0);
		dg_rec_electron.addDataSet(H_FT_e_p_f, 1);
		dg_rec_electron.addDataSet(H_FT_e_t_f, 2);
		dg_rec_electron.addDataSet(hi_rec_e_dp, 3);
		dg_rec_electron.addDataSet(fn_rec_e_dp, 3);
		dg_rec_electron.addDataSet(hi_rec_e_dtheta, 4);
		dg_rec_electron.addDataSet(fn_rec_e_dtheta, 4);
		dg_rec_electron.addDataSet(hi_rec_e_dphi, 5);
		dg_rec_electron.addDataSet(fn_rec_e_dphi, 5);
		dg_rec_electron.addDataSet(H_FT_Q2, 6);
		dg_rec_electron.addDataSet(H_FT_W, 7);
		dg_rec_electron.addDataSet(H_virphoton, 8);
		dg_rec_electron.addDataSet(H_FT_W_Q2, 9);
		dg_rec_electron.addDataSet(hi_rec_e_dvz, 11);
		dg_rec_electron.addDataSet(fn_rec_e_dvz, 11);
		dg_rec_electron.addDataSet(hi_mc_e_p_the, 12);
		dg_rec_electron.addDataSet(hi_mc_e_p_phi, 13);
		dg_rec_electron.addDataSet(hi_mc_e_the_phi, 14);

		//e resolution
		dg_rec_e_resolution = new DataGroup(4,4);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dp_p, 0);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dp_theta, 1);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dp_phi, 2);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dp_vz, 3);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dtheta_p, 4);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dtheta_theta, 5);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dtheta_phi, 6);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dtheta_vz, 7);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dphi_p, 8);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dphi_theta, 9);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dphi_phi, 10);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dphi_vz, 11);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dvz_p, 12);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dvz_theta, 13);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dvz_phi, 14);
		dg_rec_e_resolution.addDataSet(hi_rec_e_dvz_vz, 15);

		// rec fast kp (kp1)
		dg_rec_kp1 = new DataGroup(3,4);
		dg_rec_kp1.addDataSet(hi_rec_kp1_p_the, 0);
		dg_rec_kp1.addDataSet(hi_rec_kp1_p_phi, 1);
		dg_rec_kp1.addDataSet(hi_rec_kp1_the_phi, 2);
		dg_rec_kp1.addDataSet(hi_rec_kp1_dp, 3);
		dg_rec_kp1.addDataSet(fn_rec_kp1_dp, 3);
		dg_rec_kp1.addDataSet(hi_rec_kp1_dtheta, 4);
		dg_rec_kp1.addDataSet(fn_rec_kp1_dtheta, 4);
		dg_rec_kp1.addDataSet(hi_rec_kp1_dphi, 5);
		dg_rec_kp1.addDataSet(fn_rec_kp1_dphi, 5);
		dg_rec_kp1.addDataSet(hi_rec_kp1_dvx, 6);
		dg_rec_kp1.addDataSet(fn_rec_kp1_dvx, 6);
		dg_rec_kp1.addDataSet(hi_rec_kp1_dvy, 7);
		dg_rec_kp1.addDataSet(fn_rec_kp1_dvy, 7);
		dg_rec_kp1.addDataSet(hi_rec_kp1_dvz, 8);
		dg_rec_kp1.addDataSet(fn_rec_kp1_dvz, 8);
		dg_rec_kp1.addDataSet(hi_mc_kp1_p_the, 9);
		dg_rec_kp1.addDataSet(hi_mc_kp1_p_phi, 10);
		dg_rec_kp1.addDataSet(hi_mc_kp1_the_phi, 11);

		//kp1 resolution
		dg_rec_kp1_resolution = new DataGroup(4,4);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dp_p, 0);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dp_theta, 1);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dp_phi, 2);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dp_vz, 3);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dtheta_p, 4);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dtheta_theta, 5);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dtheta_phi, 6);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dtheta_vz, 7);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dphi_p, 8);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dphi_theta, 9);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dphi_phi, 10);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dphi_vz, 11);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dvz_p, 12);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dvz_theta, 13);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dvz_phi, 14);
		dg_rec_kp1_resolution.addDataSet(hi_rec_kp1_dvz_vz, 15);

		// rec slow kp (kp2);
		dg_rec_kp2 = new DataGroup(3,4);
		dg_rec_kp2.addDataSet(hi_rec_kp2_p_the, 0);
		dg_rec_kp2.addDataSet(hi_rec_kp2_p_phi, 1);
		dg_rec_kp2.addDataSet(hi_rec_kp2_the_phi, 2);
		dg_rec_kp2.addDataSet(hi_rec_kp2_dp, 3);
		dg_rec_kp2.addDataSet(fn_rec_kp2_dp, 3);
		dg_rec_kp2.addDataSet(hi_rec_kp2_dtheta, 4);
		dg_rec_kp2.addDataSet(fn_rec_kp2_dtheta, 4);
		dg_rec_kp2.addDataSet(hi_rec_kp2_dphi, 5);
		dg_rec_kp2.addDataSet(fn_rec_kp2_dphi, 5);
		dg_rec_kp2.addDataSet(hi_rec_kp2_dvx, 6);
		dg_rec_kp2.addDataSet(fn_rec_kp2_dvx, 6);
		dg_rec_kp2.addDataSet(hi_rec_kp2_dvy, 7);
		dg_rec_kp2.addDataSet(fn_rec_kp2_dvy, 7);
		dg_rec_kp2.addDataSet(hi_rec_kp2_dvz, 8);
		dg_rec_kp2.addDataSet(fn_rec_kp2_dvz, 8);
		dg_rec_kp2.addDataSet(hi_mc_kp2_p_the, 9);
		dg_rec_kp2.addDataSet(hi_mc_kp2_p_phi, 10);
		dg_rec_kp2.addDataSet(hi_mc_kp2_the_phi, 11);

		//kp2 resolution
		dg_rec_kp2_resolution = new DataGroup(4,4);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dp_p, 0);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dp_theta, 1);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dp_phi, 2);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dp_vz, 3);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dtheta_p, 4);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dtheta_theta, 5);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dtheta_phi, 6);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dtheta_vz, 7);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dphi_p, 8);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dphi_theta, 9);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dphi_phi, 10);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dphi_vz, 11);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dvz_p, 12);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dvz_theta, 13);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dvz_phi, 14);
		dg_rec_kp2_resolution.addDataSet(hi_rec_kp2_dvz_vz, 15);

		// rec km
		dg_rec_km = new DataGroup(3,4);
		dg_rec_km.addDataSet(hi_rec_km_p_the, 0);
		dg_rec_km.addDataSet(hi_rec_km_p_phi, 1);
		dg_rec_km.addDataSet(hi_rec_km_the_phi, 2);
		dg_rec_km.addDataSet(hi_rec_km_dp, 3);
		dg_rec_km.addDataSet(fn_rec_km_dp, 3);
		dg_rec_km.addDataSet(hi_rec_km_dtheta, 4);
		dg_rec_km.addDataSet(fn_rec_km_dtheta, 4);
		dg_rec_km.addDataSet(hi_rec_km_dphi, 5);
		dg_rec_km.addDataSet(fn_rec_km_dphi, 5);
		dg_rec_km.addDataSet(hi_rec_km_dvx, 6);
		dg_rec_km.addDataSet(fn_rec_km_dvx, 6);
		dg_rec_km.addDataSet(hi_rec_km_dvy, 7);
		dg_rec_km.addDataSet(fn_rec_km_dvy, 7);
		dg_rec_km.addDataSet(hi_rec_km_dvz, 8);
		dg_rec_km.addDataSet(fn_rec_km_dvz, 8);
		dg_rec_km.addDataSet(hi_mc_km_p_the, 9);
		dg_rec_km.addDataSet(hi_mc_km_p_phi, 10);
		dg_rec_km.addDataSet(hi_mc_km_the_phi, 11);

		dg_rec_km_resolution = new DataGroup(4,4);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dp_p, 0);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dp_theta, 1);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dp_phi, 2);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dp_vz, 3);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dtheta_p, 4);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dtheta_theta, 5);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dtheta_phi, 6);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dtheta_vz, 7);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dphi_p, 8);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dphi_theta, 9);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dphi_phi, 10);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dphi_vz, 11);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dvz_p, 12);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dvz_theta, 13);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dvz_phi, 14);
		dg_rec_km_resolution.addDataSet(hi_rec_km_dvz_vz, 15);

		//counter for particles
		dg_counter = new DataGroup(5, 3);
		dg_counter.addDataSet(hi_pip_counter, 0);
		dg_counter.addDataSet(hi_pim_counter, 1);
		dg_counter.addDataSet(hi_kp_counter, 2);
		dg_counter.addDataSet(hi_km_counter, 3);
		dg_counter.addDataSet(hi_prot_counter, 4);
		dg_counter.addDataSet(hi_fpip_counter, 5);
		dg_counter.addDataSet(hi_fpim_counter, 6);
		dg_counter.addDataSet(hi_fkp_counter, 7);
		dg_counter.addDataSet(hi_fkm_counter, 8);
		dg_counter.addDataSet(hi_fprot_counter, 9);
		dg_counter.addDataSet(hi_cpip_counter, 10);
		dg_counter.addDataSet(hi_cpim_counter, 11);
		dg_counter.addDataSet(hi_ckp_counter, 12);
		dg_counter.addDataSet(hi_ckm_counter, 13);
		dg_counter.addDataSet(hi_cprot_counter, 14);


		//rec xi, lambda L_lamda
		dg_rec_xi = new DataGroup(2,4);
		dg_rec_xi.addDataSet(hi_mc_ekpkpkm_mm_ekpkp, 0);
		dg_rec_xi.addDataSet(f1_mc_xi, 0);
		dg_rec_xi.addDataSet(hi_mc_ekpkpkm_mm_ekpkpkm, 1);
		dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkp, 2);
		dg_rec_xi.addDataSet(f1_xi, 2);
		dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkpkm, 3);
		//dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkpkm_withcorrn, 3);
		//dg_rec_xi.addDataSet(f1_lambda, 3);
		dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkp_withcorrn, 4);
		dg_rec_xi.addDataSet(fn_xi_with_cor, 4);
		dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkpkm_withcorrn, 5);
		//dg_rec_xi.addDataSet(f1_sigma, 1);
		dg_rec_xi.addDataSet(hi_ekpkpkm_IM_kmlambda, 6);
		dg_rec_xi.addDataSet(fn_im_kmlambda, 6);
		dg_rec_xi.addDataSet(hi_ekpkpkm_IM_kmsigma, 7);
		dg_rec_xi.addDataSet(fn_im_kmsigma, 7);
		//dg_rec_xi.addDataSet(L_sigma, 1);

		dg_req = new DataGroup(2, 1);
		dg_req.addDataSet(hi_ekpkpkm_MM_ekpkp, 0);
		dg_req.addDataSet(f1_xi, 0);
		dg_req.addDataSet(hi_ekpkpkm_MM_ekpkpkm, 1);

		dg_mm = new DataGroup(3,3);
		dg_mm.addDataSet(H_ekpkp_MM_req_pim, 0);
		dg_mm.addDataSet(H_efkpckp_MM_req_pim, 1);
		dg_mm.addDataSet(H_efkpckp_MM_req_cpim, 2);
		dg_mm.addDataSet(hi_ekpkpkm_MM_ekpkp, 3);
		dg_mm.addDataSet(H_efkpckp_MM_req_km, 4);
		dg_mm.addDataSet(H_efkpckp_MM_req_ckm, 5);
		dg_mm.addDataSet(hi_ekpkpkmprot_MM2, 6);
		
		// verticx time vs p
		dg_vtime = new DataGroup(5, 2);
		dg_vtime.addDataSet(H_pip_vt_p, 0);
		dg_vtime.addDataSet(H_pim_vt_p, 1);
		dg_vtime.addDataSet(H_kp_vt_p, 2);
		dg_vtime.addDataSet(H_km_vt_p, 3);
		dg_vtime.addDataSet(H_prot_vt_p, 4);
		dg_vtime.addDataSet(H_pipc_vt_p, 5);
		dg_vtime.addDataSet(H_pimc_vt_p, 6);
		dg_vtime.addDataSet(H_kpc_vt_p, 7);
		dg_vtime.addDataSet(H_kmc_vt_p, 8);
		dg_vtime.addDataSet(H_protc_vt_p, 9);
		
		dg_fdtof = new DataGroup(5, 2);
		dg_fdtof.addDataSet(H_FD_pos_beta_mom, 0);
		dg_fdtof.addDataSet(F_prot_beta_mom, 0);
		dg_fdtof.addDataSet(F_pip_beta_mom,0);
		dg_fdtof.addDataSet(F_kp_beta_mom,0);
		dg_fdtof.addDataSet(H_FD_pos_mass_mom, 1);
		dg_fdtof.addDataSet(H_FD_pos_mass_the, 2);
		dg_fdtof.addDataSet(H_FD_pos_mass_phi, 3);
		dg_fdtof.addDataSet(hi_fd_pos_mass, 4);
		dg_fdtof.addDataSet(H_FD_neg_beta_mom, 5);
		dg_fdtof.addDataSet(F_prot_beta_mom, 5);
		dg_fdtof.addDataSet(F_pip_beta_mom, 5);
		dg_fdtof.addDataSet(F_kp_beta_mom, 5);
		dg_fdtof.addDataSet(H_FD_neg_mass_mom, 6);
		dg_fdtof.addDataSet(H_FD_neg_mass_the, 7);
		dg_fdtof.addDataSet(H_FD_neg_mass_phi, 8);
		dg_fdtof.addDataSet(hi_fd_neg_mass, 9);

		dg_cdtof = new DataGroup(5, 2);
		dg_cdtof.addDataSet(H_CD_pos_beta_mom, 0);
		dg_cdtof.addDataSet(F_prot_beta_mom, 0);
		dg_cdtof.addDataSet(F_pip_beta_mom,0);
		dg_cdtof.addDataSet(F_kp_beta_mom,0);
		dg_cdtof.addDataSet(H_CD_pos_mass_mom, 1);
		dg_cdtof.addDataSet(H_CD_pos_mass_the, 2);
		dg_cdtof.addDataSet(H_CD_pos_mass_phi, 3);
		dg_cdtof.addDataSet(hi_cd_pos_mass, 4);
		dg_cdtof.addDataSet(H_CD_neg_beta_mom, 5);
		dg_cdtof.addDataSet(F_prot_beta_mom, 5);
		dg_cdtof.addDataSet(F_pip_beta_mom, 5);
		dg_cdtof.addDataSet(F_kp_beta_mom, 5);
		dg_cdtof.addDataSet(H_CD_neg_mass_mom, 6);
		dg_cdtof.addDataSet(H_CD_neg_mass_the, 7);
		dg_cdtof.addDataSet(H_CD_neg_mass_phi, 8);
		dg_cdtof.addDataSet(hi_cd_neg_mass, 9);

		// CD/FD hadron particles vz
		dg_vz = new DataGroup(5, 2);
		dg_vz.addDataSet(H_pip_vz, 0);
		dg_vz.addDataSet(H_pim_vz, 1);
		dg_vz.addDataSet(H_kp_vz, 2);
		dg_vz.addDataSet(H_km_vz, 3);
		dg_vz.addDataSet(H_prot_vz, 4);
		dg_vz.addDataSet(H_pipc_vz, 5);
		dg_vz.addDataSet(H_pimc_vz, 6);
		dg_vz.addDataSet(H_kpc_vz, 7);
		dg_vz.addDataSet(H_kmc_vz, 8);
		dg_vz.addDataSet(H_protc_vz, 9);

		//TOF_t, 

		dg_tof_t = new DataGroup(5, 2);
		dg_tof_t.addDataSet(H_pip_FTOF1b_t, 0);
		dg_tof_t.addDataSet(H_pim_FTOF1b_t, 1);
		dg_tof_t.addDataSet(H_kp_FTOF1b_t, 2);
		dg_tof_t.addDataSet(H_km_FTOF1b_t, 3);
		dg_tof_t.addDataSet(H_prot_FTOF1b_t, 4);
		dg_tof_t.addDataSet(H_pip_CTOF_t, 5);
		dg_tof_t.addDataSet(H_pim_CTOF_t, 6);
		dg_tof_t.addDataSet(H_kp_CTOF_t, 7);
		dg_tof_t.addDataSet(H_km_CTOF_t, 8);
		dg_tof_t.addDataSet(H_prot_CTOF_t, 9);

		//TOF_path
		dg_tof_path = new DataGroup(5, 2);
		dg_tof_path.addDataSet(H_pip_FTOF1b_path, 0);
		dg_tof_path.addDataSet(H_pim_FTOF1b_path, 1);
		dg_tof_path.addDataSet(H_kp_FTOF1b_path, 2);
		dg_tof_path.addDataSet(H_km_FTOF1b_path, 3);
		dg_tof_path.addDataSet(H_prot_FTOF1b_path, 4);
		dg_tof_path.addDataSet(H_pip_CTOF_path, 5);
		dg_tof_path.addDataSet(H_pim_CTOF_path, 6);
		dg_tof_path.addDataSet(H_kp_CTOF_path, 7);
		dg_tof_path.addDataSet(H_km_CTOF_path, 8);
		dg_tof_path.addDataSet(H_prot_CTOF_path, 9);


		
	} // end of ft_ana()

	public void fillMCPartBank(DataBank mcbank){

		mcels = new ArrayList<Particle>();
		mckps = new ArrayList<Particle>();
		mckms = new ArrayList<Particle>();

		Particle mcEl = null;
    	Particle mcKp = null;
    	Particle mcKm  = null;
    	Particle mcPro = null;
    	Particle mcPim = null;
    	LorentzVector xi		 = null;
   		LorentzVector lambda          = null;

		for (int k = 0; k < mcbank.rows(); k++) {
			if (mcbank.getInt("pid", k) == 11){
				mcEl = new Particle(
					mcbank.getInt("pid", k),
					mcbank.getFloat("px", k), 
					mcbank.getFloat("py", k), 
					mcbank.getFloat("pz", k), 
					mcbank.getFloat("vx", k), 
					mcbank.getFloat("vy", k), 
					mcbank.getFloat("vz", k));
				mcels.add(mcEl);

			}
			else if (mcbank.getInt("pid", k) == 321){
				mcKp = new Particle(
					mcbank.getInt("pid", k),
					mcbank.getFloat("px", k), 
					mcbank.getFloat("py", k), 
					mcbank.getFloat("pz", k), 
					mcbank.getFloat("vx", k), 
					mcbank.getFloat("vy", k), 
					mcbank.getFloat("vz", k));
				mckps.add(mcKp);
			}
			else if (mcbank.getInt("pid", k) == -321){
				mcKm = new Particle(
					mcbank.getInt("pid", k),
					mcbank.getFloat("px", k), 
					mcbank.getFloat("py", k), 
					mcbank.getFloat("pz", k), 
					mcbank.getFloat("vx", k), 
					mcbank.getFloat("vy", k), 
					mcbank.getFloat("vz", k));
				mckms.add(mcKm);
			}

		}
		
		if (mcels.size() == 1 && mckms.size() == 1 && mckps.size() == 2){
			dg_rec_electron.getH2F("hi_mc_e_p_the").fill(Math.toDegrees(mcels.get(0).theta()), mcels.get(0).p());
			dg_rec_electron.getH2F("hi_mc_e_p_phi").fill(Math.toDegrees(mcels.get(0).phi()), mcels.get(0).p());
			dg_rec_electron.getH2F("hi_mc_e_the_phi").fill(Math.toDegrees(mcels.get(0).phi()), Math.toDegrees(mcels.get(0).theta()));
			dg_rec_km.getH2F("hi_mc_km_p_the").fill(Math.toDegrees(mckms.get(0).theta()), mckms.get(0).p());
			dg_rec_km.getH2F("hi_mc_km_p_phi").fill(Math.toDegrees(mckms.get(0).phi()), mckms.get(0).p());
			dg_rec_km.getH2F("hi_mc_km_the_phi").fill(Math.toDegrees(mckms.get(0).phi()), Math.toDegrees(mckms.get(0).theta()));
			if(mckps.get(0).p() > mckps.get(1).p()){
				dg_rec_kp1.getH2F("hi_mc_kp1_p_the").fill(Math.toDegrees(mckps.get(0).theta()), mckps.get(0).p());
				dg_rec_kp1.getH2F("hi_mc_kp1_p_phi").fill(Math.toDegrees(mckps.get(0).phi()), mckps.get(0).p());
				dg_rec_kp1.getH2F("hi_mc_kp1_the_phi").fill(Math.toDegrees(mckps.get(0).phi()), Math.toDegrees(mckps.get(0).theta()));
				dg_rec_kp2.getH2F("hi_mc_kp2_p_the").fill(Math.toDegrees(mckps.get(1).theta()), mckps.get(1).p());
				dg_rec_kp2.getH2F("hi_mc_kp2_p_phi").fill(Math.toDegrees(mckps.get(1).phi()), mckps.get(1).p());
				dg_rec_kp2.getH2F("hi_mc_kp2_the_phi").fill(Math.toDegrees(mckps.get(1).phi()), Math.toDegrees(mckps.get(1).theta()));
			} else {
				dg_rec_kp1.getH2F("hi_mc_kp1_p_the").fill(Math.toDegrees(mckps.get(1).theta()), mckps.get(1).p());
				dg_rec_kp1.getH2F("hi_mc_kp1_p_phi").fill(Math.toDegrees(mckps.get(1).phi()), mckps.get(1).p());
				dg_rec_kp1.getH2F("hi_mc_kp1_the_phi").fill(Math.toDegrees(mckps.get(1).phi()), Math.toDegrees(mckps.get(1).theta()));
				dg_rec_kp2.getH2F("hi_mc_kp2_p_the").fill(Math.toDegrees(mckps.get(0).theta()), mckps.get(0).p());
				dg_rec_kp2.getH2F("hi_mc_kp2_p_phi").fill(Math.toDegrees(mckps.get(0).phi()), mckps.get(0).p());
				dg_rec_kp2.getH2F("hi_mc_kp2_the_phi").fill(Math.toDegrees(mckps.get(0).phi()), Math.toDegrees(mckps.get(0).theta()));				
			}

			xi = new LorentzVector(0.0, 0.0, Eb, 0.9383+Eb);
    		xi.sub(mcels.get(0).vector());
    		xi.sub(mckps.get(0).vector());
    		xi.sub(mckps.get(1).vector());

    		lambda = new LorentzVector(0.0, 0.0, Eb, 0.9383+Eb);
    		lambda.sub(mcels.get(0).vector());
    		lambda.sub(mckps.get(0).vector());
    		lambda.sub(mckps.get(1).vector());
    		lambda.sub(mckms.get(0).vector());

	    	dg_rec_xi.getH1F("hi_mc_ekpkpkm_mm_ekpkp").fill(xi.mass());
	    	dg_rec_xi.getH1F("hi_mc_ekpkpkm_mm_ekpkpkm").fill(lambda.mass()); 


		}


	}
	
	public void fillRecBank(DataBank recBank) {
		STT = recBank.getFloat("startTime", 0);
		//RFT = recBank.getFloat("RFTime", 0);
	}	
		
	
	public int makeFTElectron(DataBank bank, DataBank recFT, DataBank recftResponse) {
		int NFTElec = 0;
		ftels = new ArrayList<Particle>();
		Particle Vftel = null;
		//Mp = 0.93827f;
		e_ftCal_hitPosition = new Vector3D(-1000, -1000, -1000);
		for (int k = 0; k < bank.rows(); k++) {
			int ftbpid = recFT.getInt("pid", k);
			byte q = bank.getByte("charge", k);
			float px = bank.getFloat("px", k);
			float py = bank.getFloat("py", k);
			float pz = bank.getFloat("pz", k);
			float vx = bank.getFloat("vx", k);
			float vy = bank.getFloat("vy", k);
			float vz = bank.getFloat("vz", k);
			int partstatus = bank.getShort("status", k);
			if (partstatus<0) partstatus = -partstatus;
			int recftstatus = recFT.getShort("status", k);
			boolean inFT = (partstatus >= 1000 && partstatus < 2000 && recftstatus < 0.0);
			e_mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			e_the = (float) Math.toDegrees(Math.acos(pz / e_mom));
			if (ftbpid == 11 && q == -1 && inFT  && e_mom > 0.1 && e_mom < 4.5 && e_the < 4.5 && e_the > 2.5) { //&& e_mom > 0.5 && e_mom < 4.5
				NFTElec++; 
				Vftel = new Particle(ftbpid, px, py, pz, vx, vy, vz);
				ftels.add(Vftel);
			}

			if (NFTElec == 1) {
				found_eFT = true;
				int electron_ind = k;

				// store electron hit position on FT calorimeter
				for (int r = 0; r < recftResponse.rows(); r++) {
					if (recftResponse.getByte("detector", r) == 10 && recftResponse.getShort("pindex", r) == electron_ind) {
						e_ftCal_hitPosition.setXYZ(recftResponse.getFloat("x", r), recftResponse.getFloat("y", r), recftResponse.getFloat("z", r));
					}
				}

				e_phi = (float) Math.toDegrees(Math.atan2(py, px));
				
				//Geraints Correction factor
				double eft = Math.sqrt(e_mom * e_mom + PDGDatabase.getParticleById(11).mass()*PDGDatabase.getParticleById(11).mass());
				double eftNew = eft-0.03689+0.1412*eft-0.04316*eft*eft+0.007046*eft*eft*eft-0.0004055*eft*eft*eft*eft;
				double corrn = eftNew/eft;
				//Ve = new LorentzVector(px*corrn, py*corrn, pz*corrn, eft*corrn);
				
				//double corrn = -0.001385172 * e_mom * e_mom + 0.003978123 *  e_mom + 1.036788;
				
				//double corrn = 1.041 + 0.168 / Math.pow((10.789 - e_mom), 3.89);
				//double corrn = 1.057 + (0.0256)/Math.pow((9.95 - e_mom), 0.68);// 05-26-2020
				//double corrn = 1.097 - 0.0275 * e_mom + 0.0051 * Math.pow(e_mom, 2) - 0.0003 * Math.pow(e_mom, 3);// 05-26-2020
				//double corrn = 1.258 - 0.093 * e_mom + 0.014 * e_mom * e_mom - 0.0007 * e_mom *e_mom * e_mom;
				//double corrn = 0.982 + 0.017 * e_mom - 0.001 * Math.pow(e_mom, 2) - 0.0001 * Math.pow(e_mom, 3);
				//double corrn = 0.319 + 0.467 * e_mom -0.111 * Math.pow(e_mom, 2) + 0.012 * Math.pow(e_mom, 3) - 0.001 * Math.pow(e_mom, 4);
				//double corrn = -0.001 * e_mom * e_mom * e_mom + 0.018 * e_mom * e_mom - 0.11 * e_mom + 1.27;
				//double corrn = (double)1.045;
				//Ve = new LorentzVector(corrn*px, corrn*py, corrn*pz,  Math.sqrt(corrn * e_mom * corrn * e_mom + PDGDatabase.getParticleById(11).mass()*PDGDatabase.getParticleById(11).mass()));
				
				Ve = new LorentzVector(px, py, pz, eft);
				
				//Ve_consCorn = new LorentzVector(corrn*px, corrn*py, corrn*pz,  Math.sqrt(corrn * e_mom * corrn * e_mom + PDGDatabase.getParticleById(11).mass()*PDGDatabase.getParticleById(11).mass()));
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
	
	public void makeOthers(DataBank recbank,  DataBank recFTbank, DataBank recSCBank) {
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
			int ftbpid = recFTbank.getInt("pid", k);
			float px = recbank.getFloat("px", k);
			float py = recbank.getFloat("py", k);
			float pz = recbank.getFloat("pz", k);
			float vx = recbank.getFloat("vx", k);
			float vy = recbank.getFloat("vy", k);
			float vz = recbank.getFloat("vz", k);
			float be = recbank.getFloat("beta", k);
			float chi2pid = recFTbank.getFloat("chi2pid", k);
			float ftbbe = recFTbank.getFloat("beta", k);
			int status = recbank.getShort("status", k);
			if (status<0) status = -status;
			boolean inDC = (status >= 2000 && status < 4000);
			boolean inCD = (status >= 4000);
			float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			float the = (float) Math.toDegrees(Math.acos(pz / mom));
			float phi = (float) Math.toDegrees(Math.atan2(py, px));
			float FTBmass = mom * mom * (1 / ( ftbbe * ftbbe ) - 1);
			boolean fdkpMass2Cut = (0.2 < FTBmass && FTBmass < 0.40); // cut for K- candidates in FD 
			boolean cdkpMass2Cut = (0.15 < FTBmass && FTBmass < 0.35); //  cut for K- candidates in CD 
			//boolean kpMass2Cut = (0.180625 < FTBmass && FTBmass < 0.36);// 0.425 < kpM < 0.6
			boolean fdChi2pidCut = (Math.abs(chi2pid) < 5.0);
			boolean cdChi2pidCut = (Math.abs(chi2pid) < 5.0);
			if (ftbpid == 211 ) {
				npip++;
				// Mike's and Alan's momentum correction 
				//float pipCor = FTBmass/(FTBmass + Math.abs(chi2pid) * 0.0271f); //from TOF mass
				//float pipCor = (float) PDGDatabase.getParticleById(211).mass()/(float)(PDGDatabase.getParticleById(211).mass() + chi2pid * 0.0271f);   //from nominal mass
				pip_part_ind = k;
				//pip_mom = mom*pipCor; //with Alan's correction
				pip_mom = mom;
				pip_the = the;
				pip_phi = (float) Math.toDegrees(Math.atan2(py, px));
				pip_vx = vx;
				pip_vy = vy;
				pip_vz = vz;
				pip_ftb_beta = ftbbe;
				pip_status = status;
				
				for (int r = 0; r < recSCBank.rows(); r++) {
					if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {	
						if (recSCBank.getShort("pindex", r) == pip_part_ind && recSCBank.getByte("layer", r) == 2) {	
							nfdpip++;
							pip_FTOF_pad1b = recSCBank.getShort("component", r);
							pip_FTOF1b_t = recSCBank.getFloat("time", r);
							pip_FTOF1b_path = recSCBank.getFloat("path", r);
							float pip_beta = pip_mom / (float) Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass());
							//pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_ftb_beta * 29.98f) - STT - pip_vz/ (pip_ftb_beta * 29.98f);
							pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_beta * 29.98f) - STT - pip_vz/ (pip_beta * 29.98f);
							//pip_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							float pip_FTOF1b_TOF = (float) pip_FTOF1b_t - STT - pip_vz/ (pip_beta * 29.98f);
							if (Math.abs(pip_FTOF1b_vt) < 0.5) { //0.4
								nfdwithcutpip++;
								//Vpip = new LorentzVector(pip_px*pipCor, pip_py*pipCor, pip_pz*pipCor, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								//Vpip = new LorentzVector(pip_px, pip_py, pip_pz, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								Vpip = new Particle(ftbpid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);
								fpips.add(Vpip);
								pips.add(Vpip);
								//H_pip_vt_p.fill(pip_FTOF1b_vt, pip_mom);
							}
							H_pip_vt_p.fill(pip_FTOF1b_vt, pip_mom);
							H_pip_vz.fill(pip_vz);
							H_pip_FTOF1b_t.fill(pip_FTOF1b_TOF);
							H_pip_FTOF1b_path.fill(pip_FTOF1b_path);
							
						} // pip from FTOF panal 1b
					}
					
					if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut) {
						
						if (recSCBank.getShort("pindex", r) == pip_part_ind ) {
							ncdpip++;
							pip_CTOF_pad = recSCBank.getShort("component", r);
							pip_CTOF_t = recSCBank.getFloat("time", r);
							pip_CTOF_path = recSCBank.getFloat("path", r);
							float pipc_beta = pip_mom / (float) Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass());
							//pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_ftb_beta * 29.98f) - STT - pip_vz/ (pip_ftb_beta * 29.98f);
							pip_CTOF_vt = pip_CTOF_t - pip_CTOF_path / (pipc_beta * 29.98f) - STT - pip_vz/ (pipc_beta * 29.98f);
							//pip_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							if (Math.abs(pip_CTOF_vt) < 0.4) {
								ncdwithcutpip++;
								//Vpipc = new LorentzVector(pip_px*pipCor, pip_py*pipCor, pip_pz*pipCor, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								Vpipc = new Particle(ftbpid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);	
								pips.add(Vpipc);
								cpips.add(Vpipc);
								//H_pip_vt_p.fill(pip_FTOF1b_vt, pip_mom);
							}
							H_pipc_vt_p.fill(pip_CTOF_vt, pip_mom);
							H_pipc_vz.fill(pip_vz);
							H_pip_CTOF_t.fill(pip_CTOF_t - STT - pip_vz/ (pipc_beta * 29.98f));
							H_pip_CTOF_path.fill(pip_CTOF_path);
							
							
						}
						
						
					} //CTOF
					
					
				}			

			}
			if (ftbpid == -211 ) {
				// Mike's and Alan's momentum correction 
				//float pimCor = FTBmass/(FTBmass + Math.abs(chi2pid) * 0.0125f);// from TOF mass
				//float pimCor = (float)PDGDatabase.getParticleById(-211).mass()/(float)( PDGDatabase.getParticleById(-211).mass()+ chi2pid * 0.0125f); // from nominal mass
				pim_part_ind = k;
				pim_mom = mom;
				//pim_mom = mom*pimCor; //with Alan's correction
				pim_the = the;
				pim_phi = (float) Math.toDegrees(Math.atan2(py, px));
				pim_vx = vx;
				pim_vy = vy;
				pim_vz = vz;
				pim_ftb_beta = ftbbe;
				
				for (int r = 0; r < recSCBank.rows(); r++) {
					if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {	
						if (recSCBank.getShort("pindex", r) == pim_part_ind && recSCBank.getByte("layer", r) == 2) {
							pim_FTOF_pad1b = recSCBank.getShort("component", r);
							pim_FTOF1b_t = recSCBank.getFloat("time", r);
							pim_FTOF1b_path = recSCBank.getFloat("path", r);
							float pim_beta = pim_mom / (float) Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass() * PDGDatabase.getParticleById(-211).mass());
							//pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_ftb_beta * 29.98f) - STT - pim_vz/ (pim_ftb_beta * 29.98f);
							pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_beta * 29.98f) - STT - pim_vz/ (pim_beta * 29.98f);
							//pim_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							if (Math.abs(pim_FTOF1b_vt) < 0.5) {
								//Vpim = new LorentzVector(pim_px*pimCor, pim_py*pimCor, pim_pz*pimCor, Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass()));
								Vpim = new Particle(ftbpid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);
								pims.add(Vpim);
								fpims.add(Vpim);
								//H_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);
							}
							H_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);
							H_pim_vz.fill(pim_vz);
							H_pim_FTOF1b_t.fill(pim_FTOF1b_t - STT - pim_vz/ (pim_beta * 29.98f));
							H_pim_FTOF1b_path.fill(pim_FTOF1b_path);
						} // pim from FTOF panal 1b
					}
					
					if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut ) {
						
						if (recSCBank.getShort("pindex", r) == pim_part_ind ) {
							pim_CTOF_pad = recSCBank.getShort("component", r);
							pim_CTOF_t = recSCBank.getFloat("time", r);
							pim_CTOF_path = recSCBank.getFloat("path", r);
							float pimc_beta = pim_mom / (float) Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass());
							//pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_ftb_beta * 29.98f) - STT - pim_vz/ (pim_ftb_beta * 29.98f);
							pim_CTOF_vt = pim_CTOF_t - pim_CTOF_path / (pimc_beta * 29.98f) - STT - pim_vz/ (pimc_beta * 29.98f);
							//pim_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							if (Math.abs(pim_CTOF_vt) < 0.4) {
								//Vpimc = new LorentzVector(pim_px*pimCor, pim_py*pimCor, pim_pz*pimCor, Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass()));
								Vpimc = new Particle(ftbpid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);		
								pims.add(Vpimc);
								cpims.add(Vpimc);
								//H_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);
							}
							H_pimc_vt_p.fill(pim_CTOF_vt, pim_mom);
							H_pimc_vz.fill(pim_vz);
							H_pim_CTOF_t.fill(pim_CTOF_t - STT - pim_vz/ (pimc_beta * 29.98f));
							H_pim_CTOF_path.fill(pim_CTOF_path);
							
							
						}
						
						
					} //CTOF
					
					
				}
			}
			if (ftbpid == 321 ) {
				kp_part_ind = k;
				kp_mom = mom;
				kp_the = the;
				kp_phi = (float) Math.toDegrees(Math.atan2(py, px));
				kp_px = px;
				kp_py = py;
				kp_pz = pz;
				kp_vx = vx;
				kp_vy = vy;
				kp_vz = vz;
				kp_ftb_beta = ftbbe;
				
				for (int r = 0; r < recSCBank.rows(); r++) {
					if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut ) {	//&& fdkpMass2Cut
						if (recSCBank.getShort("pindex", r) == kp_part_ind && recSCBank.getByte("layer", r) == 2) {
							kp_FTOF_pad1b = recSCBank.getShort("component", r);
							kp_FTOF1b_t = recSCBank.getFloat("time", r);
							kp_FTOF1b_path = recSCBank.getFloat("path", r);
							float kp_beta = kp_mom / (float) Math.sqrt(kp_mom * kp_mom + PDGDatabase.getParticleById(321).mass() * PDGDatabase.getParticleById(321).mass());
							//kp_FTOF1b_vt = kp_FTOF1b_t - kp_FTOF1b_path / (k_ftb_beta * 29.98f) - STT - kp_vz/ (kp_ftb_beta * 29.98f);
							kp_FTOF1b_vt = kp_FTOF1b_t - kp_FTOF1b_path / (kp_beta * 29.98f) - STT - kp_vz/ (kp_beta * 29.98f);
							//kp_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							boolean fdkpvzCut = kp_vz > -10 && kp_vz < 2;
							if (Math.abs(kp_FTOF1b_vt) < 0.5 ) {// && kp_mom < 2.8 //fdkpvzCut && 
								Vkp = new Particle(ftbpid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);
								fkps.add(Vkp);
								kps.add(Vkp);
								//H_kp_vt_p.fill(kp_FTOF1b_vt, kp_mom);
							}
							H_kp_vt_p.fill(kp_FTOF1b_vt, kp_mom);
							H_kp_vz.fill(kp_vz);
							H_kp_FTOF1b_t.fill(kp_FTOF1b_t - STT - kp_vz/ (kp_beta * 29.98f));
							H_kp_FTOF1b_path.fill(kp_FTOF1b_path);
						} // kp from FTOF panal 1b
					}
					
					if (recSCBank.getByte("detector", r) == 4 && inCD  && cdChi2pidCut ) { //&& cdkpMass2Cut
						
						if (recSCBank.getShort("pindex", r) == kp_part_ind ) {
							kp_CTOF_pad = recSCBank.getShort("component", r);
							kp_CTOF_t = recSCBank.getFloat("time", r);
							kp_CTOF_path = recSCBank.getFloat("path", r);
							float kpc_beta = kp_mom / (float) Math.sqrt(kp_mom * kp_mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass());
							kp_CTOF_vt = kp_CTOF_t - kp_CTOF_path / (kpc_beta * 29.98f) - STT - kp_vz/ (kpc_beta * 29.98f);
							//kp_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							boolean cdkpvzCut = kp_vz > -8 && kp_vz < 1;
							if (Math.abs(kp_CTOF_vt) < 0.4) { //cdkpvzCut &&
								Vkpc = new Particle(ftbpid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);
								kps.add(Vkpc);
								ckps.add(Vkpc);
								//H_kp_vt_p.fill(kp_FTOF1b_vt, kp_mom);
							}
							H_kpc_vt_p.fill(kp_CTOF_vt, kp_mom);
							H_kpc_vz.fill(kp_vz);
							H_kp_CTOF_t.fill(kp_CTOF_t - STT - kp_vz/ (kpc_beta * 29.98f));
							H_kp_CTOF_path.fill(kp_CTOF_path);
							
							
						}
						
						
					} //CTOF
					
					
				}
				nkp++;

			}
			
			if (ftbpid == -321) {
				km_part_ind = k;
				km_mom = mom;
				km_the = the;
				km_phi = (float) Math.toDegrees(Math.atan2(py, px));
				km_px  = px;
				km_py  = py;
				km_pz  = pz;
				km_vx = vx;
				km_vy = vy;
				km_vz = vz;
				km_ftb_beta = ftbbe;
				
				for (int r = 0; r < recSCBank.rows(); r++) {
					if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {	
						if (recSCBank.getShort("pindex", r) == km_part_ind && recSCBank.getByte("layer", r) == 2) {
							km_FTOF_pad1b = recSCBank.getShort("component", r);
							km_FTOF1b_t = recSCBank.getFloat("time", r);
							km_FTOF1b_path = recSCBank.getFloat("path", r);
							float km_beta = km_mom / (float) Math.sqrt(km_mom * km_mom + PDGDatabase.getParticleById(-321).mass() * PDGDatabase.getParticleById(-321).mass());
							//km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_ftb_beta * 29.98f) - STT - km_vz/ (km_ftb_beta * 29.98f);
							km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_beta * 29.98f) - STT - km_vz/ (km_beta * 29.98f);
							//km_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							boolean fdkmvzCut = km_vz > -10 && km_vz < 2;
							float km_FTOF1b_TOF = (float) km_FTOF1b_t - STT - km_vz/ (km_beta * 29.98f);
							boolean km_FTOF1b_TOFCut = 22 < km_FTOF1b_TOF && km_FTOF1b_TOF < 28;
							if (Math.abs(km_FTOF1b_vt) < 0.5) {//&& fdkmvzCut && km_FTOF1b_TOFCut
								Vkm = new Particle(ftbpid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);
								fkms.add(Vkm);
								kms.add(Vkm);
								//H_km_vt_p.fill(km_FTOF1b_vt, km_mom);
							}
							H_km_vt_p.fill(km_FTOF1b_vt, km_mom);
							H_km_vz.fill(km_vz);
							H_km_FTOF1b_t.fill(km_FTOF1b_TOF);
							H_km_FTOF1b_path.fill(km_FTOF1b_path);
						} // km from FTOF panal 1b
					}
					
					if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut ) {
						
						if (recSCBank.getShort("pindex", r) == km_part_ind ) {
							km_CTOF_pad = recSCBank.getShort("component", r);
							km_CTOF_t = recSCBank.getFloat("time", r);
							km_CTOF_path = recSCBank.getFloat("path", r);
							float kmc_beta = km_mom / (float) Math.sqrt(km_mom * km_mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass());
							//km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_ftb_beta * 29.98f) - STT - km_vz/ (km_ftb_beta * 29.98f);
							km_CTOF_vt = km_CTOF_t - km_CTOF_path / (kmc_beta * 29.98f) - STT - km_vz/ (kmc_beta * 29.98f);
							//km_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							boolean cdkmvzCut = kp_vz > -8 && kp_vz < 1;
							if (Math.abs(km_CTOF_vt) < 0.4 ) { //&& cdkmvzCut
								Vkmc = new Particle(ftbpid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);		
								kms.add(Vkmc);
								ckms.add(Vkmc);
								//H_km_vt_p.fill(km_FTOF1b_vt, km_mom);
							}
							H_kmc_vt_p.fill(km_CTOF_vt, km_mom);
							H_kmc_vz.fill(km_vz);
							H_km_CTOF_t.fill(km_CTOF_t - STT - km_vz/ (kmc_beta * 29.98f));
							H_km_CTOF_path.fill(km_CTOF_path);
							
							
						}
						
						
					} //CTOF				
				}
			}
			
			  
			 
			
			if (ftbpid == 2212 ) {
				// Mike's and Alan's momentum correction 
				//float proCor1 = FTBmass/(FTBmass + Math.abs(chi2pid) * 0.00893f); //from TOF mass
				//float proCor2 = FTBmass/(FTBmass + Math.abs(chi2pid) * 0.004f);
				//float proCor = proCor1 * proCor2;
				
				//float proCor1 = (float)PDGDatabase.getParticleById(2212).mass()/(float)(PDGDatabase.getParticleById(2212).mass() + chi2pid * 0.00893f); // from nominal mass
				//float proCor2 = (float)PDGDatabase.getParticleById(2212).mass()/(float)(PDGDatabase.getParticleById(2212).mass() + chi2pid * 0.004f);
				//float proCor = proCor1 * proCor2;
				prot_part_ind = k;
				prot_mom = mom;
				//prot_mom = mom * proCor; //with Alan's correction
				prot_the = the;
				prot_phi = (float) Math.toDegrees(Math.atan2(py, px));
				prot_px = px;
				prot_py = py;
				prot_pz = pz;
				prot_vx = vx;
				prot_vy = vy;
				prot_vz = vz;
				prot_ftb_beta = ftbbe;
				
				for (int r = 0; r < recSCBank.rows(); r++) {
					if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {	
						if (recSCBank.getShort("pindex", r) == prot_part_ind && recSCBank.getByte("layer", r) == 2) {
							prot_FTOF_pad1b = recSCBank.getShort("component", r);
							prot_FTOF1b_t = recSCBank.getFloat("time", r);
							prot_FTOF1b_path = recSCBank.getFloat("path", r);
							float prot_beta = prot_mom / (float) Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass() * PDGDatabase.getParticleById(2212).mass());
							//prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_ftb_beta * 29.98f) - STT - prot_vz/ (prot_ftb_beta * 29.98f);
							prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_beta * 29.98f) - STT - prot_vz/ (prot_beta * 29.98f);
							//prot_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							float prot_FTOF1b_TOF = (float) prot_FTOF1b_t - STT - prot_vz/ (prot_beta * 29.98f);
							boolean prot_FTOF1b_TOFCut = 22 < prot_FTOF1b_TOF && prot_FTOF1b_TOF < 32;
							if (Math.abs(prot_FTOF1b_vt) < 0.5 ) {// && prot_FTOF1b_TOFCut
								//Vprot = new LorentzVector(prot_px * proCor, prot_py* proCor, prot_pz* proCor, Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass()));
								Vprot = new Particle(ftbpid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);
								fprots.add(Vprot);
								prots.add(Vprot);
								//H_prot_vt_p.fill(prot_FTOF1b_vt, prot_mom);
							}
							H_prot_vt_p.fill(prot_FTOF1b_vt, prot_mom);
							H_prot_vz.fill(prot_vz);
							H_prot_FTOF1b_t.fill(prot_FTOF1b_t - STT - prot_vz/ (prot_beta * 29.98f));
							H_prot_FTOF1b_path.fill(prot_FTOF1b_path);
						} // prot from FTOF panal 1b
					}
					
					if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut ) {
						
						if (recSCBank.getShort("pindex", r) == prot_part_ind ) {
							prot_CTOF_pad = recSCBank.getShort("component", r);
							prot_CTOF_t = recSCBank.getFloat("time", r);
							prot_CTOF_path = recSCBank.getFloat("path", r);
							float protc_beta = prot_mom / (float) Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass());
							//prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_ftb_beta * 29.98f) - STT - prot_vz/ (prot_ftb_beta * 29.98f);
							prot_CTOF_vt = prot_CTOF_t - prot_CTOF_path / (protc_beta * 29.98f) - STT - prot_vz/ (protc_beta * 29.98f);
							//prot_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							if (Math.abs(prot_CTOF_vt) < 0.4) {
								//Vprotc = new LorentzVector(prot_px * proCor, prot_py* proCor, prot_pz* proCor, Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass()));
								Vprotc = new Particle(ftbpid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);		
								cprots.add(Vprotc);
								prots.add(Vprotc);
								//H_prot_vt_p.fill(prot_FTOF1b_vt, prot_mom);
							}
							H_protc_vt_p.fill(prot_CTOF_vt, prot_mom);
							H_protc_vz.fill(prot_vz);
							H_prot_CTOF_t.fill(prot_CTOF_t - STT - prot_vz/ (protc_beta * 29.98f));
							H_prot_CTOF_path.fill(prot_CTOF_path);
							
							
						}
						
						
					} //CTOF
					
					
				}
			
			}
			
			
			//fillFTOF(recbank, recSCBank);
			
			
			
			if ( q > 0 ) {
				if (inDC && fdChi2pidCut) {
					hi_fd_pos_mass.fill(FTBmass);
					H_FD_pos_beta_mom.fill(mom, ftbbe);
					H_FD_pos_mass_mom.fill(mom, FTBmass);
					H_FD_pos_mass_the.fill(the, FTBmass);
					H_FD_pos_mass_phi.fill(phi, FTBmass);		
				}
				if (inCD && cdChi2pidCut) {
					hi_cd_pos_mass.fill(FTBmass);
					H_CD_pos_beta_mom.fill(mom, ftbbe);
					H_CD_pos_mass_mom.fill(mom, FTBmass);
					H_CD_pos_mass_the.fill(the, FTBmass);
					H_CD_pos_mass_phi.fill(phi, FTBmass);				
				}
					
			}
			if ( q < 0 ) {
				if (inDC && fdChi2pidCut) {
					hi_fd_neg_mass.fill(FTBmass);
					H_FD_neg_beta_mom.fill(mom, ftbbe);
					H_FD_neg_mass_mom.fill(mom, FTBmass);
					H_FD_neg_mass_the.fill(the, FTBmass);
					H_FD_neg_mass_phi.fill(phi, FTBmass);		
				}
				if (inCD && cdChi2pidCut) {
					hi_cd_neg_mass.fill(FTBmass);
					H_CD_neg_beta_mom.fill(mom, ftbbe);
					H_CD_neg_mass_mom.fill(mom, FTBmass);
					H_CD_neg_mass_the.fill(the, FTBmass);
					H_CD_neg_mass_phi.fill(phi, FTBmass);				
				}

			}
			if (q == 0 ) {
				if (inDC && fdChi2pidCut) {
					H_FD_neutral_beta_mom.fill(mom, ftbbe);
					H_FD_neutral_mass_mom.fill(mom, FTBmass);
					H_FD_neutral_mass_the.fill(the, FTBmass);
					H_FD_neg_mass_phi.fill(phi, FTBmass);		
				}
				if (inCD && cdChi2pidCut) {
					H_CD_neutral_beta_mom.fill(mom, ftbbe);
					H_CD_neutral_mass_mom.fill(mom, FTBmass);
					H_CD_neutral_mass_the.fill(the, FTBmass);
					H_CD_neutral_mass_phi.fill(phi, FTBmass);				
				}
			}
			
		} // FOR LOOp
		
//		if(npip >= 2 && nfdpip == 2) {
//		System.out.println("found total pip in this event :: " + npip +" with " + nfdpip + " pips tracks in FD " + nfdwithcutpip +" with cut & " + ncdpip + " pip tracks in CD " + ncdwithcutpip + " with cut" );
//		}

	
	} //MAKEOTHER
	
	public boolean select_ekpkp_req_pim() {
		boolean res = false;
		if( found_eFT && kps.size() >= 2 && pims.size() >= 1) {
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(kps.get(0).vector());Vmissekpkp.sub(kps.get(1).vector());
			ekpkp_MM_req_pim = (float)Vmissekpkp.mass();
			res = true;
		}
		return res;
		
	}
	
	public boolean select_efkpckp_req_pim() {
		boolean res = false;
		if( found_eFT && fkps.size() >= 1  && ckps.size() >= 1 && pims.size() >= 1) {
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			efkpckp_MM_req_pim = (float)Vmissekpkp.mass();
			res = true;
		}
		return res;
	}
	
	public boolean select_efkpckp_req_cpim() {
		boolean res = false;
		if( found_eFT && fkps.size() >= 1  && ckps.size() >= 1 && cpims.size() >= 1) {
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			efkpckp_MM_req_cpim = (float)Vmissekpkp.mass();
			res = true;
		}
		return res;
	}
	 
	public boolean select_efkpckp_req_ckm() {
		boolean res = false;
		if( found_eFT && fkps.size() >= 1  && ckps.size() >= 1 && ckms.size() >= 1) {
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			efkpckp_MM_req_ckm = (float)Vmissekpkp.mass();
			res = true;
		}
		return res;
	}
	
	public boolean select_efkpckp_req_km() {
		boolean res = false;
		if( found_eFT && fkps.size() >= 1  && ckps.size() >= 1 && kms.size() >= 1) {
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			efkpckp_MM_req_km = (float)Vmissekpkp.mass();
			res = true;
		}
		return res;
	}

	public boolean select_ekpkpkmprot(){
		boolean res = false;
		if(found_eFT && kps.size() >= 2 && kms.size() >= 1 && prots.size() >= 1){

			LorentzVector Vmissekpkpkmprot = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkmprot.add(VT);Vmissekpkpkmprot.add(VB);Vmissekpkpkmprot.sub(Ve);Vmissekpkpkmprot.sub(kps.get(0).vector());Vmissekpkpkmprot.sub(kps.get(1).vector());Vmissekpkpkmprot.sub(kms.get(0).vector());Vmissekpkpkmprot.sub(prots.get(0).vector());
			ekpkpkmprot_MM2 = (float)Vmissekpkpkmprot.mass2();
			res = true;
		}
		return res;
	}
	
	// electron kp kp and km detected 
	public boolean select_ekpkp_req_km(){
		boolean res = false;
		//ftels = new ArrayList<Particle>();
		if( found_eFT && kps.size() >= 2 && kms.size() >= 1 && mcels.size() ==1 && mckps.size() == 2 && mckms.size() == 1) {
			// labeling reconstructed kps with momentum
			if (mckps.get(0).p() > mckps.get(1).p()){
				mckp1 = mckps.get(0);
				mckp2 = mckps.get(1);
			} else {
				mckp1 = mckps.get(1);
				mckp2 = mckps.get(0);
			}

			// labeling reconstructed kps with momentum
			if (kps.get(0).p() > kps.get(1).p()){
				reckp1 = kps.get(0);
				reckp2 = kps.get(1);
			} else {
				reckp1 = kps.get(1);
				reckp2 = kps.get(0);
			}

			rec_kp1_p = (float)reckp1.p();
			rec_kp1_the = (float)Math.toDegrees(reckp1.theta());
			rec_kp1_phi = (float)Math.toDegrees(reckp1.phi());
			rec_kp1_vz = (float)reckp1.vz();

			rec_kp2_p = (float)reckp2.p();
			rec_kp2_the = (float)Math.toDegrees(reckp2.theta());
			rec_kp2_phi = (float)Math.toDegrees(reckp2.phi());
			rec_kp2_vz = (float)reckp2.vz();

			rec_km_p = (float)kms.get(0).p();
			rec_km_the = (float)Math.toDegrees(kms.get(0).theta());
			rec_km_phi = (float)Math.toDegrees(kms.get(0).phi());
			rec_km_vz = (float)kms.get(0).vz();

			// common vz for kp, kp , km and ft electron
			double avg_vz = (double)(rec_kp1_vz + rec_kp2_vz + rec_km_vz)/3;
			//mcels.get(0).vz()

			// with common vertex correction applied to the electron
			LorentzVector Ve_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);

			//LorentzVector Ve_corrected = electron_corn(Ve, e_ftCal_hitPosition, mcels.get(0).vz());
			LorentzVector Vmissekpkp_withcorn = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp_withcorn.add(VT);Vmissekpkp_withcorn.add(VB);Vmissekpkp_withcorn.sub(Ve_corrected);Vmissekpkp_withcorn.sub(reckp1.vector());Vmissekpkp_withcorn.sub(reckp2.vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			ekpkpkm_MM_ekpkp_withcorrn = (float)Vmissekpkp_withcorn.mass();

			LorentzVector Vmissekpkpkm_withcorn = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm_withcorn.add(VT);Vmissekpkpkm_withcorn.add(VB);Vmissekpkpkm_withcorn.sub(Ve_corrected);Vmissekpkpkm_withcorn.sub(reckp1.vector());Vmissekpkpkm_withcorn.sub(reckp2.vector());Vmissekpkpkm_withcorn.sub(kms.get(0).vector());
			ekpkpkm_MM_ekpkpkm_withcorrn = (float)Vmissekpkpkm_withcorn.mass();


			// without the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(reckp1.vector());Vmissekpkp.sub(reckp2.vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			ekpkpkm_MM_ekpkp = (float)Vmissekpkp.mass();



			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Ve);Vmissekpkpkm.sub(reckp1.vector());Vmissekpkpkm.sub(reckp2.vector());Vmissekpkpkm.sub(kms.get(0).vector());
			ekpkpkm_MM_ekpkpkm = (float)Vmissekpkpkm.mass();

			if (Math.abs(ekpkpkm_MM_ekpkpkm - 1.115683) <= 0.052){
				found_Lambda = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedlambda = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.115683f*1.115683f));
				LorentzVector VxifromLambda = new LorentzVector(0, 0, 0, 0);
				VxifromLambda.add(Vconstrainedlambda);VxifromLambda.add(kms.get(0).vector());
				ekpkpkm_IM_kmlambda = (float)VxifromLambda.mass();
				
			}

			if (Math.abs(ekpkpkm_MM_ekpkpkm - 1.192642) <= 0.052){
				found_Sigma = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedsigma = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.192642f*1.192642f));
				LorentzVector VxifromSigma = new LorentzVector(0, 0, 0, 0);
				VxifromSigma.add(Vconstrainedsigma);VxifromSigma.add(kms.get(0).vector());
				ekpkpkm_IM_kmsigma = (float)VxifromSigma.mass();
				
			}
			

			res = true;
		}
		return res;
		
	}

	// function to correct FT electron four vector according to common vertics 
	public LorentzVector electron_corn(Particle el_lv, Vector3D el_FtCal_hitPosition, double vz){

		Vector3D vertex = new Vector3D(el_lv.vx(), el_lv.vy(), vz);
		Vector3D line = new Vector3D(el_FtCal_hitPosition);
		line.sub(vertex);
		double thetaCorr = Math.exp(1.797 - 4.485*el_lv.e()) + Math.exp(-0.8671 - 1.078*el_lv.e());
		thetaCorr = Math.toRadians(thetaCorr);
		double phiCorr = Math.exp(4.918 - 3.828*el_lv.e()) + Math.exp(3.841 - 1.256*el_lv.e()) + Math.exp(2.874 - 0.2195*el_lv.e());
		double field = -1; //-1 for negative inbending and +1 for negative outbending. 
		phiCorr = Math.toRadians(phiCorr * field);
		double new_px = (double) el_lv.e()*Math.sin(line.theta() + thetaCorr)*Math.cos(line.phi() - phiCorr);
		double new_py = (double) el_lv.e()*Math.sin(line.theta() + thetaCorr)*Math.sin(line.phi() - phiCorr);
		double new_pz = (double) el_lv.e()*Math.cos(line.theta() + thetaCorr);

		//double dtheta = vz/(double) Math.sqrt(el_FtCal_hitPosition.x()*el_FtCal_hitPosition.x() + el_FtCal_hitPosition.y()*el_FtCal_hitPosition.y() + (el_FtCal_hitPosition.z() - vz)*(el_FtCal_hitPosition.z() - vz));
		//double theta = (double) Math.toDegrees(el_lv.theta());
		//double phi = (double) Math.toDegrees(el_lv.phi());
		//double new_px = (double) el_lv.p()*Math.sin(Math.toRadians(theta+dtheta))*Math.cos(el_lv.phi());
		//double new_py = (double) el_lv.p()*Math.sin(Math.toRadians(theta+dtheta))*Math.sin(el_lv.phi());
		//double new_pz = (double) el_lv.p()*Math.cos(Math.toRadians(theta+dtheta));
		//double new_p = (double) Math.sqrt(new_px * new_px + new_py * new_py + new_pz * new_pz);

	//This calculates a new LorentzVector for the electron assuming a straight line from the new vertex to the hit in the FTCAL.

		// double r = Math.sqrt(el_FtCal_hitPosition.x()*el_FtCal_hitPosition.x() + el_FtCal_hitPosition.y()*el_FtCal_hitPosition.y() + (el_FtCal_hitPosition.z() - vz)*(el_FtCal_hitPosition.z() - vz));
		// double px = (double) el_lv.p()*(el_FtCal_hitPosition.x()/r);
		// double py = (double) el_lv.p()*(el_FtCal_hitPosition.y()/r);
		// double pz = (double) el_lv.p()*((el_FtCal_hitPosition.z() - vz)/r);
		// LorentzVector ftCal_Ve = new LorentzVector(0, 0, 0, 0);
		// ftCal_Ve.setPxPyPzM(px, py, pz, PDGDatabase.getParticleById(11).mass()); // use this if you don't want to swim electron from FTCal hit position to common vertex in magnetic field. 

	// In order to calculate the angles at the production vertex, one needs to swim the electron back to the vertex. This is done by the following two functions for theta and phi
		
		//double theta_corn = Math.toDegrees(Math.exp(1.797 - 4.485*ftCal_Ve.e()) + Math.exp(-0.8671 - 1.078*ftCal_Ve.e()));
		// double theta_corn = Math.exp(1.797 - 4.485*ftCal_Ve.e()) + Math.exp(-0.8671 - 1.078*ftCal_Ve.e()); // small electron theta correction at production vertex in radian 
		// double theta = (double) Math.toDegrees(ftCal_Ve.theta());
		// //double phi_corn = Math.toDegrees(Math.exp(4.918 - 3.828*ftCal_Ve.e()) + Math.exp(3.841 - 1.256*ftCal_Ve.e()) + Math.exp(2.874 - 0.2195*ftCal_Ve.e()));
		// double phi_corn = Math.exp(4.918 - 3.828*ftCal_Ve.e()) + Math.exp(3.841 - 1.256*ftCal_Ve.e()) + Math.exp(2.874 - 0.2195*ftCal_Ve.e());
		// double phi = (double) Math.toDegrees(ftCal_Ve.phi());

		// //Assumes inbending data correct phi is (phi + phi_corn) for outbending it would be the (phi - phi_corn):
		// double new_px = (double) ftCal_Ve.p()*Math.sin(Math.toRadians(theta + theta_corn))*Math.cos(Math.toRadians(phi + phi_corn));
		// double new_py = (double) ftCal_Ve.p()*Math.sin(Math.toRadians(theta + theta_corn))*Math.sin(Math.toRadians(phi + phi_corn));
		// double new_pz = (double) ftCal_Ve.p()*Math.cos(Math.toRadians(theta + theta_corn));
		// LorentzVector new_Ve = new LorentzVector(0, 0, 0, 0);
		// new_Ve.setPxPyPzM(new_px, new_py, new_pz, PDGDatabase.getParticleById(11).mass());
	// theta correction 
		//el_FtCal_hitPosition.setZ(el_FtCal_hitPosition.z() - vz);
		//double theta_corn_1 = (el_lv.theta() * el_FtCal_hitPosition.z() * vz)/el_FtCal_hitPosition.r2();
		//double new_px = (double) el_lv.p()*Math.sin(el_lv.theta() + theta_corn_1 )*Math.cos(el_lv.phi());
		//double new_py = (double) el_lv.p()*Math.sin(el_lv.theta() + theta_corn_1)*Math.sin(el_lv.phi());
		//double new_pz = (double) el_lv.p()*Math.cos(el_lv.theta() + theta_corn_1);
		LorentzVector new_Ve = new LorentzVector(0, 0, 0, 0);
		new_Ve.setPxPyPzE(new_px, new_py, new_pz, el_lv.e());

		return new_Ve; 
	}
	

	public void resetCounters() {
		e_ft_part_ind = -1;
		pip_part_ind = -1;
		pim_part_ind = -1;
		kp_part_ind = -1;
		km_part_ind = -1;
		prot_part_ind = -1;
		pip_FTOF_pad1b = -1;
		pim_FTOF_pad1b = -1;
		kp_FTOF_pad1b = -1;
		km_FTOF_pad1b = -1;
		found_eFT = false;
		found_Lambda = false;
		found_Sigma = false;
		found_Cascade = false;
		found_recPip = false;
		mcels = null;
		ftels = null;
		mckps = null;
		mckp1 = null;
		mckp2 = null;
		mckms = null;
		kps = null;
		reckp1 = null;
		reckp2 = null;
		pips = null;
		kms = null;
		pims = null;
		prots = null;
		ckps = null;
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
			
	}

	public void processEvent(DataEvent event) {
		resetCounters();

		if (event.hasBank("RECFT::Event"))
			fillRecBank(event.getBank("RECFT::Event"));
		if (event.hasBank("REC::Particle") && event.hasBank("RECFT::Particle") && event.hasBank("REC::ForwardTagger")) e_ft_part_ind = makeFTElectron(event.getBank("REC::Particle"), event.getBank("RECFT::Particle"), event.getBank("REC::ForwardTagger"));
		if (event.hasBank("MC::Particle") == true) fillMCPartBank(event.getBank("MC::Particle"));
		if(e_ft_part_ind > -1 && found_eFT) {		
			if (event.hasBank("REC::Particle") && event.hasBank("RECFT::Particle") && event.hasBank("REC::Scintillator")) makeOthers(event.getBank("REC::Particle"), event.getBank("RECFT::Particle"), event.getBank("REC::Scintillator"));
			
			//if (event.hasBank("REC::Scintillator")) fillFTOF(event.getBank("REC::Particle"), event.getBank("REC::Scintillator"));
			//System.out.println("found pip :: " + pips.size() + " pips tracks in FD & " + cpips.size() + " pip tracks in CD.");
			FillHists();
			
		} // e_ft_part_ind > -1
		
		
	} //processEvent
	
	public void analyze() {
		fitxi(dg_rec_xi.getH1F("hi_ekpkpkm_MM_ekpkp_withcorrn"), dg_rec_xi.getF1D("fn_xi_with_cor"));
    	fitxi(dg_rec_xi.getH1F("hi_ekpkpkm_MM_ekpkp"), dg_rec_xi.getF1D("f1_xi"));
    	fitxi(dg_rec_xi.getH1F("hi_mc_ekpkpkm_mm_ekpkp"), dg_rec_xi.getF1D("f1_mc_xi"));
    	fitxi(dg_rec_xi.getH1F("hi_ekpkpkm_IM_kmlambda"), dg_rec_xi.getF1D("fn_im_kmlambda"));
    	fitxi(dg_rec_xi.getH1F("hi_ekpkpkm_IM_kmsigma"), dg_rec_xi.getF1D("fn_im_kmsigma"));

    	fit_dp(dg_rec_electron.getH1F("hi_rec_e_dp"), dg_rec_electron.getF1D("fn_rec_e_dp"));
    	fit_del(dg_rec_electron.getH1F("hi_rec_e_dtheta"), dg_rec_electron.getF1D("fn_rec_e_dtheta"));
    	fit_del(dg_rec_electron.getH1F("hi_rec_e_dphi"), dg_rec_electron.getF1D("fn_rec_e_dphi"));
    	fit_del(dg_rec_electron.getH1F("hi_rec_e_dvz"), dg_rec_electron.getF1D("fn_rec_e_dvz"));

    	fit_dp(dg_rec_kp1.getH1F("hi_rec_kp1_dp"), dg_rec_kp1.getF1D("fn_rec_kp1_dp"));
    	fit_del(dg_rec_kp1.getH1F("hi_rec_kp1_dtheta"), dg_rec_kp1.getF1D("fn_rec_kp1_dtheta"));
    	fit_del(dg_rec_kp1.getH1F("hi_rec_kp1_dphi"), dg_rec_kp1.getF1D("fn_rec_kp1_dphi"));
    	fit_del(dg_rec_kp1.getH1F("hi_rec_kp1_dvx"), dg_rec_kp1.getF1D("fn_rec_kp1_dvx"));
    	fit_del(dg_rec_kp1.getH1F("hi_rec_kp1_dvy"), dg_rec_kp1.getF1D("fn_rec_kp1_dvy"));
    	fit_del(dg_rec_kp1.getH1F("hi_rec_kp1_dvz"), dg_rec_kp1.getF1D("fn_rec_kp1_dvz"));

    	fit_dp(dg_rec_kp2.getH1F("hi_rec_kp2_dp"), dg_rec_kp2.getF1D("fn_rec_kp2_dp"));
    	fit_del(dg_rec_kp2.getH1F("hi_rec_kp2_dtheta"), dg_rec_kp2.getF1D("fn_rec_kp2_dtheta"));
    	fit_del(dg_rec_kp2.getH1F("hi_rec_kp2_dphi"), dg_rec_kp2.getF1D("fn_rec_kp2_dphi"));
    	fit_del(dg_rec_kp2.getH1F("hi_rec_kp2_dvx"), dg_rec_kp2.getF1D("fn_rec_kp2_dvx"));
    	fit_del(dg_rec_kp2.getH1F("hi_rec_kp2_dvy"), dg_rec_kp2.getF1D("fn_rec_kp2_dvy"));
    	fit_del(dg_rec_kp2.getH1F("hi_rec_kp2_dvz"), dg_rec_kp2.getF1D("fn_rec_kp2_dvz"));

    	fit_dp(dg_rec_km.getH1F("hi_rec_km_dp"), dg_rec_km.getF1D("fn_rec_km_dp"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dtheta"), dg_rec_km.getF1D("fn_rec_km_dtheta"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dphi"), dg_rec_km.getF1D("fn_rec_km_dphi"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dvx"), dg_rec_km.getF1D("fn_rec_km_dvx"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dvy"), dg_rec_km.getF1D("fn_rec_km_dvy"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dvz"), dg_rec_km.getF1D("fn_rec_km_dvz"));
    	

    	//fitlambda(dg_rec_xi.getH1F("hi_ekpkpkm_MM_ekpkpkm"), dg_rec_xi.getF1D("f1_lambda"));
    	//fitlambda(dg_rec_xi.getH1F("hi_ekpkpkm_MM_ekpkpkm"), dg_rec_xi.getF1D("f1_sigma"));
    	//fitvz(dg_rec_kp1.getH1F("hi_rec_kp1_dvz"), dg_rec_kp1.getF1D("f_gaus"));
    	//fitvz(dg_rec_kp2.getH1F("hi_rec_kp2_dvz"), dg_rec_kp2.getF1D("f_gaus"));
    	//fitvz(dg_rec_km.getH1F("hi_rec_km_dvz"), dg_rec_km.getF1D("f_gaus"));
    //	fitxi(dg_mm.getH1F("H_"))


	}


	public void fitxi(H1F hixi, F1D f1xi){

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
        double rmax = mean + 3.0 * Math.abs(sigma);
        double rmin = mean - 3.0 * Math.abs(sigma);
        f1xi.setRange(rmin, rmax);
        DataFitter.fit(f1xi, hixi, "Q"); //No options uses error for sigma 
        hixi.setFunction(null);
        mean = f1xi.getParameter(1);
        sigma = f1xi.getParameter(2);
        rmax = mean + 2 * Math.abs(sigma);
        rmin = mean - 2 * Math.abs(sigma);
        f1xi.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1xi, hixi, "Q"); //No options uses error for sigma 
        hixi.setFunction(null);

	}

	public void fitlambda(H1F hiw, F1D f1w) {

		// get histogram maximum in the rane 1.08 -1.14
		int i1 = hiw.getXaxis().getBin(1.0);
		int i2 = hiw.getXaxis().getBin(1.14);
		int imaxx = hiw.getXaxis().getBin(1.05);
		int himaxy = hiw.getBinContent(imaxx);
		double hiMax = 0;
		int imax = i1;
		for(int i = i1; i <= i2; i2++){
			if (hiMax < hiw.getBinContent(i)){
				imax = i;
				hiMax = hiw.getBinContent(i);
			}
		}

     //   double mean = hiw.getDataX(imax);//hiw.getDataX(hiw.getMaximumBin());
     //   double amp = hiMax;//hiw.getBinContent(hiw.getMaximumBin());
        double mean = hiw.getDataX(imaxx);//hiw.getDataX(hiw.getMaximumBin());
        double amp = himaxy;//hiw.getBinContent(hiw.getMaximumBin());
        double sigma = 0.01;
        f1w.setParameter(0, amp);
        f1w.setParameter(1, mean);
        f1w.setParameter(2, sigma);
        double rmax = mean + 0.5 * Math.abs(sigma);
        double rmin = mean - 0.5 * Math.abs(sigma);
        //f1w.setRange(rmin, rmax);
        f1w.setRange(1.04, 1.14);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
        mean = f1w.getParameter(1);
        sigma = f1w.getParameter(2);
        rmax = mean + 1.0 * Math.abs(sigma);
        rmin = mean - 1.0 * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
    }

    public void fitsigma(H1F hiw, F1D f1w) {

		// get histogram maximum in the rane 1.17 -1.24
		int i1 = hiw.getXaxis().getBin(1.15);
		int i2 = hiw.getXaxis().getBin(1.3);
		int imaxx = hiw.getXaxis().getBin(1.19);
		int himaxy = hiw.getBinContent(imaxx);
		double hiMax = 0;
		int imax = i1;
		for(int i = i1; i <= i2; i2++){
			if (hiMax < hiw.getBinContent(i)){
				imax = i;
				hiMax = hiw.getBinContent(i);
			}
		}

     //   double mean = hiw.getDataX(imax);//hiw.getDataX(hiw.getMaximumBin());
     //   double amp = hiMax;//hiw.getBinContent(hiw.getMaximumBin());
     	double mean = hiw.getDataX(imaxx);//hiw.getDataX(hiw.getMaximumBin());
        double amp = himaxy;
        double sigma = 0.01;
        f1w.setParameter(0, amp);
        f1w.setParameter(1, mean);
        f1w.setParameter(2, sigma);
        double rmax = mean + 0.5 * Math.abs(sigma);
        double rmin = mean - 0.5 * Math.abs(sigma);
        //f1w.setRange(rmin, rmax);
        f1w.setRange(1.16, 1.24);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
        mean = f1w.getParameter(1);
        sigma = f1w.getParameter(2);
        rmax = mean + 1.0 * Math.abs(sigma);
        rmin = mean - 1.0 * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
    }

    public void fit_del(H1F hiw, F1D f1w){

    	double mean = hiw.getDataX(hiw.getMaximumBin());
        double amp = hiw.getBinContent(hiw.getMaximumBin());
        double sigma = 0.4;
        f1w.setParameter(0, amp);
        f1w.setParameter(1, mean);
        f1w.setParameter(2, sigma);
        double rmax = mean + 2.0 * Math.abs(sigma);
        double rmin = mean - 2.0 * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
        mean = f1w.getParameter(1);
        sigma = f1w.getParameter(2);
        rmax = mean + 3.0 * Math.abs(sigma);
        rmin = mean - 3.0 * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
    }
	
	public void fit_dp(H1F hiw, F1D f1w){

    	double mean = hiw.getDataX(hiw.getMaximumBin());
        double amp = hiw.getBinContent(hiw.getMaximumBin());
        double sigma = 0.04;
        f1w.setParameter(0, amp);
        f1w.setParameter(1, mean);
        f1w.setParameter(2, sigma);
        double rmax = mean + 2.0 * Math.abs(sigma);
        double rmin = mean - 2.0 * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
        mean = f1w.getParameter(1);
        sigma = f1w.getParameter(2);
        rmax = mean + 3.0 * Math.abs(sigma);
        rmin = mean - 3.0 * Math.abs(sigma);
        f1w.setRange(rmin, rmax);
//        System.out.println(mean + " " + sigma + " " + rmin + " " + rmax);
        DataFitter.fit(f1w, hiw, "Q"); //No options uses error for sigma 
        hiw.setFunction(null);
    }
	
	
	public void FillHists() {
		
		if(found_eFT){		
//			H_FT_e_t_f.fill(e_phi, e_the);
//			H_FT_e_p_f.fill(e_phi, e_mom);
//			H_FT_e_p_the.fill(e_the, e_mom);
			H_FT_W_Q2.fill(e_W, e_Q2);
			H_FT_W.fill(e_W);
			H_FT_Q2.fill(e_Q2);	
			H_virphoton.fill(e_virphoton);	
		}

		// fill particle counters with one electron detected in FT

		
		if(select_ekpkp_req_pim()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			H_ekpkp_MM_req_pim.fill(ekpkp_MM_req_pim);
		}
		if(select_efkpckp_req_pim()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			H_efkpckp_MM_req_pim.fill(efkpckp_MM_req_pim);
		}
		
		if(select_efkpckp_req_cpim()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			H_efkpckp_MM_req_cpim.fill(efkpckp_MM_req_cpim);
		}
		
		if(select_efkpckp_req_km()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			H_efkpckp_MM_req_km.fill(efkpckp_MM_req_km);
		}
		
		if(select_efkpckp_req_ckm()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			H_efkpckp_MM_req_ckm.fill(efkpckp_MM_req_ckm);
		}

		if(select_ekpkpkmprot()){
			hi_ekpkpkmprot_MM2.fill(ekpkpkmprot_MM2);
		}
		
		if(select_ekpkp_req_km()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			hi_ekpkpkm_MM_ekpkp.fill(ekpkpkm_MM_ekpkp);
			hi_ekpkpkm_MM_ekpkpkm.fill(ekpkpkm_MM_ekpkpkm);

			//with electron common vertex correction
			hi_ekpkpkm_MM_ekpkp_withcorrn.fill(ekpkpkm_MM_ekpkp_withcorrn);
			hi_ekpkpkm_MM_ekpkpkm_withcorrn.fill(ekpkpkm_MM_ekpkpkm_withcorrn);

			//IM of (kmlambd) with 2sigma cut on MM(ekpkpkm)
			if(found_Lambda){
				hi_ekpkpkm_IM_kmlambda.fill(ekpkpkm_IM_kmlambda);
			}
			
			if(found_Sigma){
				hi_ekpkpkm_IM_kmsigma.fill(ekpkpkm_IM_kmsigma);
			}
			
			//hi_ekpkpkm_IM_kmsigma.fill(ekpkpkm_IM_kmsigma);
			//IM of (kmsima) with 3sigma cut on MM(ekpkpkm)

			//average vz for kp, kp and km
			float avg_vz = (rec_kp1_vz + rec_kp2_vz + rec_km_vz)/3;

			H_FT_e_t_f.fill(e_phi, e_the);
			H_FT_e_p_f.fill(e_phi, e_mom);
			H_FT_e_p_the.fill(e_the, e_mom);
			hi_rec_e_dp.fill((e_mom - (float)mcels.get(0).p())/e_mom);
			hi_rec_e_dtheta.fill(e_the - (float)Math.toDegrees(mcels.get(0).theta()));
			hi_rec_e_dphi.fill(e_phi - (float)Math.toDegrees(mcels.get(0).phi()));
			//dg_rec_electron.getH1F("hi_rec_e_dvx").fill(avg_vx - (float)mcels.get(0).vx());
			//dg_rec_electron.getH1F("hi_rec_e_dvy").fill(avg_vy - (float)mcels.get(0).vy());
			dg_rec_electron.getH1F("hi_rec_e_dvz").fill(avg_vz - (float)mcels.get(0).vz());
			dg_rec_e_resolution.getH2F("hi_rec_e_dp_p").fill(e_mom, (e_mom - (float)mcels.get(0).p())/e_mom);
			dg_rec_e_resolution.getH2F("hi_rec_e_dp_theta").fill(e_the, (e_mom - (float)mcels.get(0).p())/e_mom);
			dg_rec_e_resolution.getH2F("hi_rec_e_dp_phi").fill(e_phi, (e_mom - (float)mcels.get(0).p())/e_mom);
			dg_rec_e_resolution.getH2F("hi_rec_e_dp_vz").fill(avg_vz, (e_mom - (float)mcels.get(0).p())/e_mom);
			dg_rec_e_resolution.getH2F("hi_rec_e_dtheta_p").fill(e_mom, e_the - (float)Math.toDegrees(mcels.get(0).theta()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dtheta_theta").fill(e_the, e_the - (float)Math.toDegrees(mcels.get(0).theta()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dtheta_phi").fill(e_phi, e_the - (float)Math.toDegrees(mcels.get(0).theta()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dtheta_vz").fill(avg_vz, e_the - (float)Math.toDegrees(mcels.get(0).theta()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dphi_p").fill(e_mom, e_phi - (float)Math.toDegrees(mcels.get(0).phi()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dphi_theta").fill(e_the, e_phi - (float)Math.toDegrees(mcels.get(0).phi()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dphi_phi").fill(e_phi, e_phi - (float)Math.toDegrees(mcels.get(0).phi()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dphi_vz").fill(avg_vz, e_phi - (float)Math.toDegrees(mcels.get(0).phi()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dvz_p").fill(e_mom, avg_vz - (float)mcels.get(0).vz());
			dg_rec_e_resolution.getH2F("hi_rec_e_dvz_theta").fill(e_the, avg_vz - (float)mcels.get(0).vz());
			dg_rec_e_resolution.getH2F("hi_rec_e_dvz_phi").fill(e_phi, avg_vz - (float)mcels.get(0).vz());
			dg_rec_e_resolution.getH2F("hi_rec_e_dvz_vz").fill(avg_vz, avg_vz - (float)mcels.get(0).vz());


			
			//e resolution
			

			// kp1 resolution
			hi_rec_kp1_p_the.fill(rec_kp1_the, rec_kp1_p);
			hi_rec_kp1_p_phi.fill(rec_kp1_phi, rec_kp1_p);
			hi_rec_kp1_the_phi.fill(rec_kp1_phi, rec_kp1_the);
			hi_rec_kp1_dp.fill((rec_kp1_p - (float)mckp1.p())/rec_kp1_p);
			hi_rec_kp1_dtheta.fill(rec_kp1_the - (float)Math.toDegrees(mckp1.theta()));
			hi_rec_kp1_dphi.fill(rec_kp1_phi - (float)Math.toDegrees(mckp1.phi()));
			hi_rec_kp1_dvx.fill((float)reckp1.vx() - (float)mckp1.vx());
			hi_rec_kp1_dvy.fill((float)reckp1.vy() - (float)mckp1.vy());
			hi_rec_kp1_dvz.fill(avg_vz - (float)mckp1.vz());

			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dp_p").fill(rec_kp1_p, (rec_kp1_p - (float)mckp1.p())/rec_kp1_p);
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dp_theta").fill(rec_kp1_the, (rec_kp1_p - (float)mckp1.p())/rec_kp1_p);
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dp_phi").fill(rec_kp1_phi, (rec_kp1_p - (float)mckp1.p())/rec_kp1_p);
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dp_vz").fill(avg_vz, (rec_kp1_p - (float)mckp1.p())/rec_kp1_p);
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dtheta_p").fill(rec_kp1_p, rec_kp1_the - (float)Math.toDegrees(mckp1.theta()));
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dtheta_theta").fill(rec_kp1_the, rec_kp1_the - (float)Math.toDegrees(mckp1.theta()));
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dtheta_phi").fill(rec_kp1_phi, rec_kp1_the - (float)Math.toDegrees(mckp1.theta()));
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dtheta_vz").fill(avg_vz, rec_kp1_the - (float)Math.toDegrees(mckp1.theta()));
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dphi_p").fill(rec_kp1_p, rec_kp1_phi - (float)Math.toDegrees(mckp1.phi()));
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dphi_theta").fill(rec_kp1_the, rec_kp1_phi - (float)Math.toDegrees(mckp1.phi()));
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dphi_phi").fill(rec_kp1_phi, rec_kp1_phi - (float)Math.toDegrees(mckp1.phi()));
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dphi_vz").fill(avg_vz, rec_kp1_phi - (float)Math.toDegrees(mckp1.phi()));
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dvz_p").fill(rec_kp1_p, avg_vz - (float)mckp1.vz());
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dvz_theta").fill(rec_kp1_the, avg_vz - (float)mckp1.vz());
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dvz_phi").fill(rec_kp1_phi, avg_vz - (float)mckp1.vz());
			dg_rec_kp1_resolution.getH2F("hi_rec_kp1_dvz_vz").fill(avg_vz, avg_vz - (float)mckp1.vz());

			//kp2 resolution
			hi_rec_kp2_p_the.fill(rec_kp2_the, rec_kp2_p);
			hi_rec_kp2_p_phi.fill(rec_kp2_phi, rec_kp2_p);
			hi_rec_kp2_the_phi.fill(rec_kp2_phi, rec_kp2_the);
			hi_rec_kp2_dp.fill((rec_kp2_p - (float)mckp2.p())/rec_kp2_p);
			hi_rec_kp2_dtheta.fill(rec_kp2_the - (float)Math.toDegrees(mckp2.theta()));
			hi_rec_kp2_dphi.fill(rec_kp2_phi - (float)Math.toDegrees(mckp2.phi()));
			hi_rec_kp2_dvx.fill((float)reckp2.vx() - (float)mckp2.vx());
			hi_rec_kp2_dvy.fill((float)reckp2.vy() - (float)mckp2.vy());
			hi_rec_kp2_dvz.fill(avg_vz - (float)mckp2.vz());

			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dp_p").fill(rec_kp2_p, (rec_kp2_p - (float)mckp2.p())/rec_kp2_p);
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dp_theta").fill(rec_kp2_the, (rec_kp2_p - (float)mckp2.p())/rec_kp2_p);
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dp_phi").fill(rec_kp2_phi, (rec_kp2_p - (float)mckp2.p())/rec_kp2_p);
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dp_vz").fill(avg_vz, (rec_kp2_p - (float)mckp2.p())/rec_kp2_p);
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dtheta_p").fill(rec_kp2_p, rec_kp2_the - (float)Math.toDegrees(mckp2.theta()));
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dtheta_theta").fill(rec_kp2_the, rec_kp2_the - (float)Math.toDegrees(mckp2.theta()));
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dtheta_phi").fill(rec_kp2_phi, rec_kp2_the - (float)Math.toDegrees(mckp2.theta()));
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dtheta_vz").fill(avg_vz, rec_kp2_the - (float)Math.toDegrees(mckp2.theta()));
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dphi_p").fill(rec_kp2_p, rec_kp2_phi - (float)Math.toDegrees(mckp2.phi()));
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dphi_theta").fill(rec_kp2_the, rec_kp2_phi - (float)Math.toDegrees(mckp2.phi()));
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dphi_phi").fill(rec_kp2_phi, rec_kp2_phi - (float)Math.toDegrees(mckp2.phi()));
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dphi_vz").fill(avg_vz, rec_kp2_phi - (float)Math.toDegrees(mckp2.phi()));
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dvz_p").fill(rec_kp2_p, avg_vz - (float)mckp2.vz());
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dvz_theta").fill(rec_kp2_the, avg_vz - (float)mckp2.vz());
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dvz_phi").fill(rec_kp2_phi, avg_vz - (float)mckp2.vz());
			dg_rec_kp2_resolution.getH2F("hi_rec_kp2_dvz_vz").fill(avg_vz, avg_vz - (float)mckp2.vz());

			//km resolution
			hi_rec_km_p_the.fill(rec_km_the, rec_km_p);
			hi_rec_km_p_phi.fill(rec_km_phi, rec_km_p);
			hi_rec_km_the_phi.fill(rec_km_phi, rec_km_the);
			hi_rec_km_dp.fill((rec_km_p - (float)mckms.get(0).p())/rec_km_p);
			hi_rec_km_dtheta.fill(rec_km_the - (float)Math.toDegrees(mckms.get(0).theta()));
			hi_rec_km_dphi.fill(rec_km_phi - (float)Math.toDegrees(mckms.get(0).phi()));
			hi_rec_km_dvx.fill((float)kms.get(0).vx() - (float)mckms.get(0).vx());
			hi_rec_km_dvy.fill((float)kms.get(0).vy() - (float)mckms.get(0).vy());
			hi_rec_km_dvz.fill(avg_vz - (float)mckms.get(0).vz());

			dg_rec_km_resolution.getH2F("hi_rec_km_dp_p").fill(rec_km_p, (rec_km_p - (float)mckms.get(0).p())/rec_km_p);
			dg_rec_km_resolution.getH2F("hi_rec_km_dp_theta").fill(rec_km_the, (rec_km_p - (float)mckms.get(0).p())/rec_km_p);
			dg_rec_km_resolution.getH2F("hi_rec_km_dp_phi").fill(rec_km_phi, (rec_km_p - (float)mckms.get(0).p())/rec_km_p);
			dg_rec_km_resolution.getH2F("hi_rec_km_dp_vz").fill(avg_vz, (rec_km_p - (float)mckms.get(0).p())/rec_km_p);
			dg_rec_km_resolution.getH2F("hi_rec_km_dtheta_p").fill(rec_km_p, rec_km_the - (float)Math.toDegrees(mckms.get(0).theta()));
			dg_rec_km_resolution.getH2F("hi_rec_km_dtheta_theta").fill(rec_km_the, rec_km_the - (float)Math.toDegrees(mckms.get(0).theta()));
			dg_rec_km_resolution.getH2F("hi_rec_km_dtheta_phi").fill(rec_km_phi, rec_km_the - (float)Math.toDegrees(mckms.get(0).theta()));
			dg_rec_km_resolution.getH2F("hi_rec_km_dtheta_vz").fill(avg_vz, rec_km_the - (float)Math.toDegrees(mckms.get(0).theta()));
			dg_rec_km_resolution.getH2F("hi_rec_km_dphi_p").fill(rec_km_p, rec_km_phi - (float)Math.toDegrees(mckms.get(0).phi()));
			dg_rec_km_resolution.getH2F("hi_rec_km_dphi_theta").fill(rec_km_the, rec_km_phi - (float)Math.toDegrees(mckms.get(0).phi()));
			dg_rec_km_resolution.getH2F("hi_rec_km_dphi_phi").fill(rec_km_phi, rec_km_phi - (float)Math.toDegrees(mckms.get(0).phi()));
			dg_rec_km_resolution.getH2F("hi_rec_km_dphi_vz").fill(avg_vz, rec_km_phi - (float)Math.toDegrees(mckms.get(0).phi()));
			//dg_rec_km_resolution.getH2F("hi_rec_km_dphi_vz").fill((float)mckms.get(0).vz(), rec_km_phi - (float)Math.toDegrees(mckms.get(0).phi()));
			dg_rec_km_resolution.getH2F("hi_rec_km_dvz_p").fill(rec_km_p, avg_vz - (float)mckms.get(0).vz());
			dg_rec_km_resolution.getH2F("hi_rec_km_dvz_theta").fill(rec_km_the, avg_vz - (float)mckms.get(0).vz());
			dg_rec_km_resolution.getH2F("hi_rec_km_dvz_phi").fill(rec_km_phi, avg_vz - (float)mckms.get(0).vz());
			dg_rec_km_resolution.getH2F("hi_rec_km_dvz_vz").fill(avg_vz, avg_vz - (float)mckms.get(0).vz());

			if(pips.size() > 0) hi_pip_counter.fill(pips.size());
			if(pims.size() > 0) hi_pim_counter.fill(pims.size());
			if(kps.size() > 0) hi_kp_counter.fill(kps.size());
			if(kms.size() > 0) hi_km_counter.fill(kms.size());
			if(prots.size() > 0) hi_prot_counter.fill(prots.size());
			if(fpips.size() > 0) hi_fpip_counter.fill(fpips.size());
			if(fpims.size() > 0) hi_fpim_counter.fill(fpims.size());
			if(fkps.size() > 0) hi_fkp_counter.fill(fkps.size());
			if(fkms.size() > 0) hi_fkm_counter.fill(fkms.size());
			if(fprots.size() > 0) hi_fprot_counter.fill(fprots.size());
			if(cpips.size() > 0) hi_cpip_counter.fill(cpips.size());
			if(cpims.size() > 0) hi_cpim_counter.fill(cpims.size());
			if(ckps.size() > 0) hi_ckp_counter.fill(ckps.size());
			if(ckms.size() > 0) hi_ckm_counter.fill(ckms.size());
			if(cprots.size() > 0) hi_cprot_counter.fill(cprots.size());

		}
		
	}
	
	
	public void plot() {

		//calculation from FD-TOF 

		myCanvas = new EmbeddedCanvasTabbed("REC-electron","electron-res","REC-kp1(fast)","kp1(fast-kp)-res","REC-kp2(slow)","kp2(slow-kp)-res", "REC-km","km-res","REC-Xi","MM-spectra", "FD-TOF", "CD-TOF", "VTime", "particle-Vz","TOF-t","TOF-path", "Counter", "Request");
		
		//reconstructed electron

		myCanvas.getCanvas("REC-electron").divide(3, 5);
		myCanvas.getCanvas("REC-electron").setGridX(false);
		myCanvas.getCanvas("REC-electron").setGridY(false);
		myCanvas.getCanvas("REC-electron").setAxisFontSize(18);
		myCanvas.getCanvas("REC-electron").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-electron").draw(dg_rec_electron);
		myCanvas.getCanvas("REC-electron").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(9).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(12).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(13).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(14).getAxisZ().setLog(true);

		myCanvas.getCanvas("electron-res").divide(4,4);
		myCanvas.getCanvas("electron-res").setGridX(false);
		myCanvas.getCanvas("electron-res").setGridY(false);
		myCanvas.getCanvas("electron-res").setAxisFontSize(18);
		myCanvas.getCanvas("electron-res").setAxisTitleSize(24);
		myCanvas.getCanvas("electron-res").draw(dg_rec_e_resolution);
		myCanvas.getCanvas("electron-res").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(4).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(9).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(10).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(11).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(12).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(13).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(14).getAxisZ().setLog(true);
		myCanvas.getCanvas("electron-res").getPad(15).getAxisZ().setLog(true);

		//reconstructed kp (kp1;fast k+ and slow k+; kp2)

		myCanvas.getCanvas("REC-kp1(fast)").divide(3,4);
		myCanvas.getCanvas("REC-kp1(fast)").setGridX(false);
		myCanvas.getCanvas("REC-kp1(fast)").setGridY(false);
		myCanvas.getCanvas("REC-kp1(fast)").setAxisFontSize(18);
		myCanvas.getCanvas("REC-kp1(fast)").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-kp1(fast)").draw(dg_rec_kp1);
		myCanvas.getCanvas("REC-kp1(fast)").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp1(fast)").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp1(fast)").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp1(fast)").getPad(9).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp1(fast)").getPad(10).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp1(fast)").getPad(11).getAxisZ().setLog(true);

		// fast kaon resolution
		myCanvas.getCanvas("kp1(fast-kp)-res").divide(4,4);
		myCanvas.getCanvas("kp1(fast-kp)-res").setGridX(false);
		myCanvas.getCanvas("kp1(fast-kp)-res").setGridY(false);
		myCanvas.getCanvas("kp1(fast-kp)-res").setAxisFontSize(18);
		myCanvas.getCanvas("kp1(fast-kp)-res").setAxisTitleSize(24);
		myCanvas.getCanvas("kp1(fast-kp)-res").draw(dg_rec_kp1_resolution);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(4).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(9).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(10).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(11).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(12).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(13).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(14).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp1(fast-kp)-res").getPad(15).getAxisZ().setLog(true);

		// slow kaon
		myCanvas.getCanvas("REC-kp2(slow)").divide(3,4);
		myCanvas.getCanvas("REC-kp2(slow)").setGridX(false);
		myCanvas.getCanvas("REC-kp2(slow)").setGridY(false);
		myCanvas.getCanvas("REC-kp2(slow)").setAxisFontSize(18);
		myCanvas.getCanvas("REC-kp2(slow)").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-kp2(slow)").draw(dg_rec_kp2);
		myCanvas.getCanvas("REC-kp2(slow)").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp2(slow)").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp2(slow)").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp2(slow)").getPad(9).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp2(slow)").getPad(10).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp2(slow)").getPad(11).getAxisZ().setLog(true);

		//slow kaon resolution
		myCanvas.getCanvas("kp2(slow-kp)-res").divide(4,4);
		myCanvas.getCanvas("kp2(slow-kp)-res").setGridX(false);
		myCanvas.getCanvas("kp2(slow-kp)-res").setGridY(false);
		myCanvas.getCanvas("kp2(slow-kp)-res").setAxisFontSize(18);
		myCanvas.getCanvas("kp2(slow-kp)-res").setAxisTitleSize(24);
		myCanvas.getCanvas("kp2(slow-kp)-res").draw(dg_rec_kp2_resolution);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(4).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(9).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(10).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(11).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(12).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(13).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(14).getAxisZ().setLog(true);
		myCanvas.getCanvas("kp2(slow-kp)-res").getPad(15).getAxisZ().setLog(true);

		//reconstructed km
		myCanvas.getCanvas("REC-km").divide(3,4);
		myCanvas.getCanvas("REC-km").setGridX(false);
		myCanvas.getCanvas("REC-km").setGridY(false);
		myCanvas.getCanvas("REC-km").setAxisFontSize(18);
		myCanvas.getCanvas("REC-km").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-km").draw(dg_rec_km);
		myCanvas.getCanvas("REC-km").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-km").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-km").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-km").getPad(9).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-km").getPad(10).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-km").getPad(11).getAxisZ().setLog(true);

		myCanvas.getCanvas("km-res").divide(4,4);
		myCanvas.getCanvas("km-res").setGridX(false);
		myCanvas.getCanvas("km-res").setGridY(false);
		myCanvas.getCanvas("km-res").setAxisFontSize(18);
		myCanvas.getCanvas("km-res").setAxisTitleSize(24);
		myCanvas.getCanvas("km-res").draw(dg_rec_km_resolution);
		myCanvas.getCanvas("km-res").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(4).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(9).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(10).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(11).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(12).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(13).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(14).getAxisZ().setLog(true);
		myCanvas.getCanvas("km-res").getPad(15).getAxisZ().setLog(true);


		myCanvas.getCanvas("REC-Xi").divide(2,4);
		myCanvas.getCanvas("REC-Xi").setGridX(false);
		myCanvas.getCanvas("REC-Xi").setGridY(false);
		myCanvas.getCanvas("REC-Xi").setAxisFontSize(18);
		myCanvas.getCanvas("REC-Xi").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-Xi").draw(dg_rec_xi);
		//myCanvas.getCanvas("REC-Xi").getPad(1).draw(L_lamda);
		//myCanvas.getCanvas("REC-Xi").getPad(1).draw(L_sigma);

		myCanvas.getCanvas("MM-spectra").divide(3,3);
		myCanvas.getCanvas("MM-spectra").setGridX(false);
		myCanvas.getCanvas("MM-spectra").setGridY(false);
		myCanvas.getCanvas("MM-spectra").setAxisFontSize(18);
		myCanvas.getCanvas("MM-spectra").setAxisTitleSize(24);
		myCanvas.getCanvas("MM-spectra").draw(dg_mm);


		myCanvas.getCanvas("FD-TOF").divide(5, 2);
		myCanvas.getCanvas("FD-TOF").setGridX(false);
		myCanvas.getCanvas("FD-TOF").setGridY(false);
		myCanvas.getCanvas("FD-TOF").setAxisFontSize(18);
		myCanvas.getCanvas("FD-TOF").setAxisTitleSize(24);
		myCanvas.getCanvas("FD-TOF").draw(dg_fdtof);
		myCanvas.getCanvas("FD-TOF").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(4).getAxisY().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(9).getAxisY().setLog(true);
		//calculations from CD-TOF
		myCanvas.getCanvas("CD-TOF").divide(5, 2);
		myCanvas.getCanvas("CD-TOF").setGridX(false);
		myCanvas.getCanvas("CD-TOF").setGridY(false);
		myCanvas.getCanvas("CD-TOF").setAxisFontSize(18);
		myCanvas.getCanvas("CD-TOF").setAxisTitleSize(24);
		myCanvas.getCanvas("CD-TOF").draw(dg_cdtof);
		myCanvas.getCanvas("CD-TOF").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(4).getAxisY().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(9).getAxisY().setLog(true);

		// vertex time hadrons
		myCanvas.getCanvas("VTime").divide(5, 2);
		myCanvas.getCanvas("VTime").setGridX(false);
		myCanvas.getCanvas("VTime").setGridY(false);
		myCanvas.getCanvas("VTime").setAxisFontSize(18);
		myCanvas.getCanvas("VTime").setAxisTitleSize(24);
		myCanvas.getCanvas("VTime").draw(dg_vtime);
		myCanvas.getCanvas("VTime").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(4).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("VTime").getPad(9).getAxisZ().setLog(true);

		// particle verticx

		myCanvas.getCanvas("particle-Vz").divide(5, 2);
		myCanvas.getCanvas("particle-Vz").setGridX(false);
		myCanvas.getCanvas("particle-Vz").setGridY(false);
		myCanvas.getCanvas("particle-Vz").setAxisFontSize(18);
		myCanvas.getCanvas("particle-Vz").setAxisTitleSize(24);
		myCanvas.getCanvas("particle-Vz").draw(dg_vz);

		myCanvas.getCanvas("TOF-t").divide(5, 2);
		myCanvas.getCanvas("TOF-t").setGridX(false);
		myCanvas.getCanvas("TOF-t").setGridY(false);
		myCanvas.getCanvas("TOF-t").setAxisFontSize(18);
		myCanvas.getCanvas("TOF-t").setAxisTitleSize(24);
		myCanvas.getCanvas("TOF-t").draw(dg_tof_t);

		myCanvas.getCanvas("TOF-path").divide(5, 2);
		myCanvas.getCanvas("TOF-path").setGridX(false);
		myCanvas.getCanvas("TOF-path").setGridY(false);
		myCanvas.getCanvas("TOF-path").setAxisFontSize(18);
		myCanvas.getCanvas("TOF-path").setAxisTitleSize(24);
		myCanvas.getCanvas("TOF-path").draw(dg_tof_path);

		
		myCanvas.getCanvas("Counter").divide(5, 3);
		myCanvas.getCanvas("Counter").setGridX(false);
		myCanvas.getCanvas("Counter").setGridY(false);
		myCanvas.getCanvas("Counter").setAxisFontSize(18);
		myCanvas.getCanvas("Counter").setAxisTitleSize(24);
		myCanvas.getCanvas("Counter").draw(dg_counter);


		myCanvas.getCanvas("Request").divide(2, 1);
		myCanvas.getCanvas("Request").setGridX(false);
		myCanvas.getCanvas("Request").setGridY(false);
		myCanvas.getCanvas("Request").setAxisFontSize(18);
		myCanvas.getCanvas("Request").setAxisTitleSize(24);
		myCanvas.getCanvas("Request").draw(dg_req);
		
				
	}
	
	public void save() {
		
		//TDirectory dirout = new TDirectory();
		//dirout.mkdir("");
		//dirout.cd("/FTElec/");
		//dirout.addDataSet(H_FT_e_t_f, H_FT_e_p_f, H_FT_W_Q2);
		myCanvas.getCanvas("REC-electron").save("ft_electron.png");
		myCanvas.getCanvas("electron-res").save("ft_electron_resolution.png");
		myCanvas.getCanvas("REC-kp1(fast)").save("fast_kaon_plus.png");
		myCanvas.getCanvas("REC-kp2(slow)").save("slow_kaon_plus.png");
		myCanvas.getCanvas("kp1(fast-kp)-res").save("fast_kaon_plus_resolution.png");
		myCanvas.getCanvas("kp2(slow-kp)-res").save("slow_kaon_plus_resolution.png");
		myCanvas.getCanvas("REC-km").save("kaon_minus.png");
		myCanvas.getCanvas("km-res").save("kaon_minus_resolution.png");
		myCanvas.getCanvas("REC-Xi").save("cascade.png");
		myCanvas.getCanvas("MM-spectra").save("mm_ekpkpkm.png");
		myCanvas.getCanvas("CD-TOF").save("cd_beta_p_overview.png");
		myCanvas.getCanvas("FD-TOF").save("fd_beta_p_overview.png");
		myCanvas.getCanvas("VTime").save("particle_dtime_overview.png");
		myCanvas.getCanvas("particle-Vz").save("particle_vz_overview.png");
		myCanvas.getCanvas("TOF-t").save("tof_time_overview.png");
		myCanvas.getCanvas("TOF-path").save("tof_path_overview.png");
		myCanvas.getCanvas("Counter").save("particle_multiplicity.png");
		myCanvas.getCanvas("Request").save("reconstructed_cascade_1820.png");

	}

	public void showplots() {

		JFrame frame = new JFrame("SIMULATION");
		frame.setSize(1600, 1000);
		frame.add(myCanvas);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	}
	
		
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	public static void main(String[] args) {
		//System.setProperty("java.awt.headless", "true");
		GStyle.setPalette("kRainBow");
		GStyle.getH1FAttributes().setOptStat("1110");
        //GStyle.getFunctionAttributes().setOptStat("1100");

		int count = 0;
		int maxevents = 10000000;
		ekpkpreqkmavgvz ana = new ekpkpreqkmavgvz();
		System.out.println(String.format(">>> files from list %s >>>", args[0]));
		String filelist = "list_of_files.txt";
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
            	DataEvent event = reader.getNextEvent();
                ana.processEvent(event);
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
        
		 		
		System.out.println("Total events : " + count);
		ana.analyze();
		ana.plot();
		ana.showplots();
		//ana.save();
		System.out.println("Good Bye !!!!!");
	}
		
	
	
}
