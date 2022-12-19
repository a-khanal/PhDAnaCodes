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

public class efdky_data_v1 {
	//public int NFTElec;
	
	public float  Eb, Mp;
	public float rec_EB, rec_beam_px, rec_beam_py, rec_beam_pz;
	//public float EB, Eb, Mp;
	public float STT, RFT, FTSTT, vt;
	
	public Particle Vprot, Vpip, Vpim, Vkp, Vkm;
	public Particle Vprotc, Vpipc, Vpimc, Vkpc, Vkmc;
	public Particle Vftel_corrected;
	public LorentzVector Ve_consCorn;
	public LorentzVector VB, VT, Ve, VGS, VhadronSystm, Vpim_correct, Vkp_correct;

	public Vector3D e_ftCal_hitPosition;
	
	public boolean found_eFD, found_eFT, found_Lambda, found_Sigma, found_Cascade, found_recPip;
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
	public List<Particle> fdels;
	public List<Particle> ftels, pips, pims, kps, kms, prots;
	public List<Particle> mcels, mckps, mckms;

    public Particle mckp1;
    public Particle mckp2;
	// for e kp kp km detected
	public Particle reckp1;
    public Particle reckp2;
	public float rec_kp_p, rec_kp_the, rec_kp_phi, rec_kp_vz, rec_prot_p, rec_prot_the, rec_prot_phi, rec_prot_vz;
	public float rec_km_p, rec_km_the, rec_km_phi, rec_km_vz;

	public float ekpprot_mm_ekp, ekpprot_mm2_ekpprot, ekpprot_mm2_ekpprot_mm_ekp;
	public float ekpkpkm_MM_ekpkp, ekpkpkm_MM_ekpkpkm, ekpkpkm_IM_kmlambda, ekpkpkm_IM_kmsigma;
	public float ekpkpkm_MM_ekpkp_nocorr, ekpkpkm_MM_ekpkpkm_nocorr;
	public float ekpkpkmprot_MM2;

	
	public F1D fn_rec_e_dp, fn_rec_e_dtheta, fn_rec_e_dphi, fn_rec_e_dvx, fn_rec_e_dvy, fn_rec_e_dvz;
	public F1D fn_rec_kp_dp, fn_rec_kp_dtheta, fn_rec_kp_dphi, fn_rec_kp_dvx, fn_rec_kp_dvy, fn_rec_kp_dvz;
	public F1D fn_rec_prot_dp, fn_rec_prot_dtheta, fn_rec_prot_dphi, fn_rec_prot_dvx, fn_rec_prot_dvy, fn_rec_prot_dvz;
	public F1D fn_rec_km_dp, fn_rec_km_dtheta, fn_rec_km_dphi, fn_rec_km_dvx, fn_rec_km_dvy, fn_rec_km_dvz;

	public H2F H_FT_e_beta_mom;
	public H2F H_FT_e_t_f, H_FT_e_p_f, H_FT_e_p_the;
	public H1F hi_rec_e_dp, hi_rec_e_dtheta, hi_rec_e_dphi;
	public H2F hi_rec_e_dp_p, hi_rec_e_dp_theta, hi_rec_e_dp_phi, hi_rec_e_dp_vz, hi_rec_e_dtheta_p, hi_rec_e_dtheta_theta, hi_rec_e_dtheta_phi, hi_rec_e_dtheta_vz, hi_rec_e_dphi_p, hi_rec_e_dphi_theta, hi_rec_e_dphi_phi, hi_rec_e_dphi_vz, hi_rec_e_dvz_p, hi_rec_e_dvz_theta, hi_rec_e_dvz_phi, hi_rec_e_dvz_vz;
	public H2F H_FT_W_Q2, H_FT_e_xB_Q2;

	// momentum, theta and phi resolution for kp1, kp2 and km
	public H1F hi_rec_kp_dp, hi_rec_kp_dtheta, hi_rec_kp_dphi;
	public H2F hi_rec_kp_dp_p, hi_rec_kp_dp_theta, hi_rec_kp_dp_phi, hi_rec_kp_dp_vz, hi_rec_kp_dtheta_p, hi_rec_kp_dtheta_theta, hi_rec_kp_dtheta_phi, hi_rec_kp_dtheta_vz, hi_rec_kp_dphi_p, hi_rec_kp_dphi_theta, hi_rec_kp_dphi_phi, hi_rec_kp_dphi_vz, hi_rec_kp_dvz_p, hi_rec_kp_dvz_theta, hi_rec_kp_dvz_phi, hi_rec_kp_dvz_vz;
	public H1F hi_rec_prot_dp, hi_rec_prot_dtheta, hi_rec_prot_dphi;
	public H2F hi_rec_prot_dp_p, hi_rec_prot_dp_theta, hi_rec_prot_dp_phi, hi_rec_prot_dp_vz, hi_rec_prot_dtheta_p, hi_rec_prot_dtheta_theta, hi_rec_prot_dtheta_phi, hi_rec_prot_dtheta_vz, hi_rec_prot_dphi_p, hi_rec_prot_dphi_theta, hi_rec_prot_dphi_phi, hi_rec_prot_dphi_vz, hi_rec_prot_dvz_p, hi_rec_prot_dvz_theta, hi_rec_prot_dvz_phi, hi_rec_prot_dvz_vz;
	public H1F hi_rec_km_dp, hi_rec_km_dtheta, hi_rec_km_dphi;
	public H2F hi_rec_km_dp_p, hi_rec_km_dp_theta, hi_rec_km_dp_phi, hi_rec_km_dp_vz, hi_rec_km_dtheta_p, hi_rec_km_dtheta_theta, hi_rec_km_dtheta_phi, hi_rec_km_dtheta_vz, hi_rec_km_dphi_p, hi_rec_km_dphi_theta, hi_rec_km_dphi_phi, hi_rec_km_dphi_vz, hi_rec_km_dvz_p, hi_rec_km_dvz_theta, hi_rec_km_dvz_phi, hi_rec_km_dvz_vz;

	// scatter plotes for hadrons 
	public H2F hi_rec_kp_p_the, hi_rec_prot_p_the, hi_rec_km_p_the, hi_mc_e_p_the, hi_mc_kp_p_the, hi_mc_prot_p_the, hi_mc_km_p_the;
	public H2F hi_rec_kp_p_phi, hi_rec_prot_p_phi, hi_rec_km_p_phi, hi_mc_e_p_phi, hi_mc_kp_p_phi, hi_mc_prot_p_phi, hi_mc_km_p_phi;
	public H2F hi_rec_kp_the_phi, hi_rec_prot_the_phi, hi_rec_km_the_phi, hi_mc_e_the_phi, hi_mc_kp_the_phi, hi_mc_prot_the_phi, hi_mc_km_the_phi;
	public H2F hi_rec_kp_theta_vz, hi_rec_prot_theta_vz;

	public H1F hi_rec_kp_vz, hi_rec_prot_vz;
	//vertex resolution for hadrons 
	public H1F hi_rec_e_dvx, hi_rec_kp_dvx, hi_rec_prot_dvx, hi_rec_km_dvx;
	public H1F hi_rec_e_dvy, hi_rec_kp_dvy, hi_rec_prot_dvy, hi_rec_km_dvy;
	public H1F hi_rec_e_dvz, hi_rec_kp_dvz, hi_rec_prot_dvz, hi_rec_km_dvz;

	//mass spectrum
	public H1F hi_ekpprot_mm2_ekpprot, hi_ekpprot_mm_ekp;

	public H2F hi_ekpprot_mm2_ekpprot_mm_ekp;

	//mass spectrum
	public H1F hi_ekpkpkm_MM_ekpkp, hi_ekpkpkm_MM_ekpkpkm, hi_ekpkpkm_IM_kmlambda, hi_ekpkpkm_IM_kmsigma;
	public H1F hi_ekpkpkm_MM_ekpkp_nocorr, hi_ekpkpkm_MM_ekpkpkm_nocorr
	public F1D f1_xi, fn_xi_no_corr, f1_mc_xi, f1_lambda, f1_sigma, f_gaus, fn_im_kmlambda, fn_im_kmsigma;
	public H1F hi_mc_ekpkpkm_mm_ekpkp, hi_mc_ekpkpkm_mm_ekpkpkm;


	public H1F hi_ekpkpkmprot_MM2;

	public H1F hi_pip_counter, hi_pim_counter, hi_kp_counter, hi_km_counter, hi_prot_counter, hi_fpip_counter, hi_fpim_counter, hi_fkp_counter, hi_fkm_counter, hi_fprot_counter, hi_cpip_counter, hi_cpim_counter, hi_ckp_counter, hi_ckm_counter, hi_cprot_counter;

	public EmbeddedCanvasTabbed myCanvas;
	public DataGroup dg_kinematics, dg_rec_electron, dg_rec_kp, dg_rec_prot, dg_rec_km, dg_rec_y, dg_rec_p, dg_rec_pim, dg_vtime, dg_fdtof, dg_cdtof, dg_vz, dg_tof_t, dg_tof_path, dg_counter;
	public DataGroup dg_rec_e_resolution, dg_rec_kp_resolution, dg_rec_prot_resolution, dg_rec_km_resolution;
	public DataGroup dg_mm2;
	public DataGroup dg_req;

	
	
	public efdky_data_v1() {

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

		//reconstructed
		H_FT_e_t_f = new H2F("H_FT_e_t_f", "H_FT_e_t_f", 100, -180, 180, 100, 0, 45);
		H_FT_e_t_f.setTitle("electron #theta vs #phi");
		H_FT_e_t_f.setTitleX("#phi (^o)");
		H_FT_e_t_f.setTitleY("#theta (^o)");
		
		H_FT_e_p_the = new H2F("H_FT_e_p_the", "H_FT_e_p_the", 100, 0, 45, 100, 0, 12);
		H_FT_e_p_the.setTitle("electron p vs #theta (^o)");
		H_FT_e_p_the.setTitleX("#theta (^o)");
		H_FT_e_p_the.setTitleY("p (GeV)");

		H_FT_e_p_f = new H2F("H_FT_e_p_f", "H_FT_e_p_f", 100, -180, 180, 100, 0, 12);
		H_FT_e_p_f.setTitle("electron p vs #phi");
		H_FT_e_p_f.setTitleX("#phi (^o)");
		H_FT_e_p_f.setTitleY("p (GeV)");

		
		
		H_FT_W_Q2 = new H2F("H_FT_W_Q2", "H_FT_W_Q2", 100, 0, 5, 100, 0, 12);
		H_FT_W_Q2.setTitle("FT Q^2 vs W");
		H_FT_W_Q2.setTitleX("W ( GeV)");
		H_FT_W_Q2.setTitleY("Q^2 (GeV^2)");
		
		
		H_FT_W = new H1F("H_FT_W", "H_FT_W", 100, 0, 5);
		H_FT_W.setTitle("electron W");
		H_FT_W.setTitleX("W (GeV)");
		H_FT_W.setTitleY("count");
		H_FT_W.setFillColor(LIGHTGREEN);


		H_FT_Q2 = new H1F("H_FT_Q2", "H_FT_Q2", 100, 0, 12);
		H_FT_Q2.setFillColor(LIGHTGREEN);
		H_FT_Q2.setTitleX("Q^2 (GeV^2)");

		H_virphoton = new H1F("H_virphoton", "H_virphoton", 100, 0, 12);
		H_virphoton.setFillColor(LIGHTGREEN);
		H_virphoton.setTitleX("E_#gamma (GeV)");


		// reconstructed
		hi_rec_kp_p_the = new H2F("hi_rec_kp_p_the", "hi_rec_kp_p_the", 100, 0, 100, 100, 0, 8);
		hi_rec_kp_p_the.setTitleX("#theta (^o)");
		hi_rec_kp_p_the.setTitleY("p (GeV)");
		hi_rec_kp_p_phi = new H2F("hi_rec_kp_p_phi", "hi_rec_kp_p_phi", 100, -180, 180, 100, 0, 8);
		hi_rec_kp_p_phi.setTitleX("#phi (^o)");
		hi_rec_kp_p_phi.setTitleY("p (GeV)");
		hi_rec_kp_the_phi = new H2F("hi_rec_kp_the_phi", "hi_rec_kp_the_phi", 100, -180, 180, 100, 0, 100);
		hi_rec_kp_the_phi.setTitleX("#phi (^o)");
		hi_rec_kp_the_phi.setTitleY("#theta (^o)");
		hi_rec_kp_theta_vz = new H2F("hi_rec_kp_theta_vz", "hi_rec_kp_theta_vz", 40, -40, 40, 100, 0, 100);
		hi_rec_kp_vz = new H1F("hi_rec_kp_vz", "hi_rec_kp_vz", 40, -40, 40);





		// reconstructed
		hi_rec_prot_p_the = new H2F("hi_rec_prot_p_the", "hi_rec_prot_p_the", 100, 0, 100, 100, 0, 8);
		hi_rec_prot_p_the.setTitleX("#theta (^o)");
		hi_rec_prot_p_the.setTitleY("p (GeV)");
		hi_rec_prot_p_phi = new H2F("hi_rec_prot_p_phi", "hi_rec_prot_p_phi", 100, -180, 180, 100, 0, 8);
		hi_rec_prot_p_phi.setTitleX("#phi (^o)");
		hi_rec_prot_p_phi.setTitleY("p (GeV)");
		hi_rec_prot_the_phi = new H2F("hi_rec_prot_the_phi", "hi_rec_prot_the_phi", 100, -180, 180, 100, 0, 100);
		hi_rec_prot_the_phi.setTitleX("#phi (^o)");
		hi_rec_prot_the_phi.setTitleY("#theta (^o)");
		hi_rec_prot_theta_vz = new H2F("hi_rec_prot_theta_vz", "hi_rec_prot_theta_vz", 40, -40, 40, 100, 0, 100);
		hi_rec_prot_vz = new H1F("hi_rec_prot_vz", "hi_rec_prot_vz", 40, -40, 40);


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
		


		// mass histos
		hi_ekpprot_mm2_ekpprot = new H1F("hi_ekpprot_mm2_ekpprot", "hi_ekpprot_mm2_ekpprot", 70, -0.2, 0.5);

		hi_ekpprot_mm_ekp = new H1F("hi_ekpprot_mm_ekp", "hi_ekpprot_mm_ekp", 50, 0.9, 1.4);

		hi_ekpprot_mm2_ekpprot_mm_ekp = new H2F("ekpprot_mm2_ekpprot_mm_ekp", "ekpprot_mm2_ekpprot_mm_ekp", 50, 0.9, 1.4, 70, -0.2, 0.5);

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
    	hi_ekpkpkm_MM_ekpkp_nocorr = new H1F("hi_ekpkpkm_MM_ekpkp_nocorr", "hi_ekpkpkm_MM_ekpkp_nocorr", 50, 1.6, 2.1);
		hi_ekpkpkm_MM_ekpkp_nocorr.setTitle("MM");
		hi_ekpkpkm_MM_ekpkp_nocorr.setTitleX("MM(eK^+K^+) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkp_nocorr.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkp_nocorr.setFillColor(LIGHTGREEN);

		fn_xi_no_corr = new F1D("fn_xi_no_corr", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		fn_xi_no_corr.setParameter(0, 0);
    	fn_xi_no_corr.setParameter(1, 1);
    	fn_xi_no_corr.setParameter(2, 0.2);
    	fn_xi_no_corr.setLineWidth(2);
    	fn_xi_no_corr.setLineColor(2);
    	fn_xi_no_corr.setOptStat("1111");

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
		hi_mc_ekpkpkm_mm_ekpkpkm.setTitleY("Events/[10 MeV/c^2]");
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
		hi_ekpkpkm_MM_ekpkpkm.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm.setFillColor(LIGHTGREEN);

		// for one with common vertex correction for electron
		hi_ekpkpkm_MM_ekpkpkm_nocorr = new H1F("hi_ekpkpkm_MM_ekpkpkm_nocorr", "hi_ekpkpkm_MM_ekpkpkm_nocorr", 50, 0.9, 1.4);
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitle("MM");
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitleX("MM(eK^+K^+K^-) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setFillColor(LIGHTGREEN);

		hi_ekpkpkm_IM_kmlambda = new H1F("hi_ekpkpkm_IM_kmlambda", "hi_ekpkpkm_IM_kmlambda", 50, 1.6, 2.1);
		hi_ekpkpkm_IM_kmlambda.setTitle("M");
		hi_ekpkpkm_IM_kmlambda.setTitleX("M(K^-#Lambda) [GeV/c^2]");
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
		hi_ekpkpkm_IM_kmsigma.setTitleX("M(K^-#Sigma) [GeV/c^2]");
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
		

		dg_kinematics = new DataGroup(2,2);
		dg_kinematics.addDataSet(H_FT_Q2, 0);
		dg_kinematics.addDataSet(H_FT_W, 1);
		dg_kinematics.addDataSet(H_FT_W_Q2, 2);

		// rec electron
		dg_rec_electron = new DataGroup(3, 2);
		dg_rec_electron.addDataSet(H_FT_e_p_the, 0);
		dg_rec_electron.addDataSet(H_FT_e_p_f, 1);
		dg_rec_electron.addDataSet(H_FT_e_t_f, 2);
		dg_rec_electron.addDataSet(H_virphoton, 3);

		// rec fast kp (kp1)
		dg_rec_kp = new DataGroup(3,2);
		dg_rec_kp.addDataSet(hi_rec_kp_p_the, 0);
		dg_rec_kp.addDataSet(hi_rec_kp_p_phi, 1);
		dg_rec_kp.addDataSet(hi_rec_kp_the_phi, 2);
		dg_rec_kp.addDataSet(hi_rec_kp_vz, 3);
		dg_rec_kp.addDataSet(hi_rec_kp_theta_vz, 4);

		// rec slow kp (kp2);
		dg_rec_prot = new DataGroup(3,2);
		dg_rec_prot.addDataSet(hi_rec_prot_p_the, 0);
		dg_rec_prot.addDataSet(hi_rec_prot_p_phi, 1);
		dg_rec_prot.addDataSet(hi_rec_prot_the_phi, 2);
		dg_rec_prot.addDataSet(hi_rec_prot_vz, 3);
		dg_rec_prot.addDataSet(hi_rec_prot_theta_vz, 4);

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


		//rec Y, lambda/sigma
		dg_rec_y = new DataGroup(2,2);
		dg_rec_y.addDataSet(hi_ekpprot_mm2_ekpprot, 0);
		dg_rec_y.addDataSet(hi_ekpprot_mm_ekp, 1);
		dg_rec_y.addDataSet(hi_ekpprot_mm2_ekpprot_mm_ekp, 2);
		//dg_rec_y.addDataSet(hi_ekpprot_mm_ekp_withcut, 3);
		//dg_rec_y.addDataSet(hi_ekpprot_mm_ekp_constraint, 4);

		/*dg_mm2 = new DataGroup(3,3);
		dg_mm2.addDataSet(hi_ekpprot_mm2_ekpprot_e_p, 0);
		dg_mm2.addDataSet(hi_ekpprot_mm2_ekpprot_e_theta, 1);
		dg_mm2.addDataSet(hi_ekpprot_mm2_ekpprot_e_phi, 2);
		dg_mm2.addDataSet(hi_ekpprot_mm2_ekpprot_kp_p, 0);
		dg_mm2.addDataSet(hi_ekpprot_mm2_ekpprot_kp_theta, 0);
		dg_mm2.addDataSet(hi_ekpprot_mm2_ekpprot_kp_phi, 0);*/

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
				dg_rec_kp.getH2F("hi_mc_kp_p_the").fill(Math.toDegrees(mckps.get(0).theta()), mckps.get(0).p());
				dg_rec_kp.getH2F("hi_mc_kp_p_phi").fill(Math.toDegrees(mckps.get(0).phi()), mckps.get(0).p());
				dg_rec_kp.getH2F("hi_mc_kp_the_phi").fill(Math.toDegrees(mckps.get(0).phi()), Math.toDegrees(mckps.get(0).theta()));
				dg_rec_prot.getH2F("hi_mc_prot_p_the").fill(Math.toDegrees(mckps.get(1).theta()), mckps.get(1).p());
				dg_rec_prot.getH2F("hi_mc_prot_p_phi").fill(Math.toDegrees(mckps.get(1).phi()), mckps.get(1).p());
				dg_rec_prot.getH2F("hi_mc_prot_the_phi").fill(Math.toDegrees(mckps.get(1).phi()), Math.toDegrees(mckps.get(1).theta()));
			} else {
				dg_rec_kp.getH2F("hi_mc_kp_p_the").fill(Math.toDegrees(mckps.get(1).theta()), mckps.get(1).p());
				dg_rec_kp.getH2F("hi_mc_kp_p_phi").fill(Math.toDegrees(mckps.get(1).phi()), mckps.get(1).p());
				dg_rec_kp.getH2F("hi_mc_kp_the_phi").fill(Math.toDegrees(mckps.get(1).phi()), Math.toDegrees(mckps.get(1).theta()));
				dg_rec_prot.getH2F("hi_mc_prot_p_the").fill(Math.toDegrees(mckps.get(0).theta()), mckps.get(0).p());
				dg_rec_prot.getH2F("hi_mc_prot_p_phi").fill(Math.toDegrees(mckps.get(0).phi()), mckps.get(0).p());
				dg_rec_prot.getH2F("hi_mc_prot_the_phi").fill(Math.toDegrees(mckps.get(0).phi()), Math.toDegrees(mckps.get(0).theta()));				
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

	    	dg_rec_y.getH1F("hi_mc_ekpkpkm_mm_ekpkp").fill(xi.mass());
	    	dg_rec_y.getH1F("hi_mc_ekpkpkm_mm_ekpkpkm").fill(lambda.mass()); 


		}


	}
	
	public void fillRecBank(DataBank recBank) {
		STT = recBank.getFloat("startTime", 0);
		//RFT = recBank.getFloat("RFTime", 0);
	}

	/*public int makeFTElectron(DataBank bank, DataBank recFT, DataBank recftResponse) {
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
			if (ftbpid == 11 && q == -1 && inFT  && e_mom > 0 && e_the < 4.5 && e_the > 2.5) { //&& e_mom > 0.5 && e_mom < 4.5
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
				//VGS.sub(Vftel.vector());
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
	}*/	

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
			e_mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			e_the = (float) Math.toDegrees(Math.acos(pz / e_mom));
			if (pid == 11 && inDC && bank.getShort("status", k) < 0.0 && q == -1){
				NFDElec++;
				Vfdel = new Particle(pid, px, py, pz, vx, vy, vz);
				fdels.add(Vfdel);
			}
			if (NFDElec == 1) {
				found_eFD = true;
				double eft = Math.sqrt(e_mom * e_mom + PDGDatabase.getParticleById(11).mass()*PDGDatabase.getParticleById(11).mass());
				e_phi = (float) Math.toDegrees(Math.atan2(py, px));
				Ve = new LorentzVector(px, py, pz, eft);
				VGS = new LorentzVector(0, 0, 0, 0);
				VGS.add(VB);
				VGS.sub(Ve);
				//VGS.sub(Vftel.vector());
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
			float mass = mom * mom * (1 / ( beta * beta ) - 1);
			//float FTBmass = mom * mom * (1 / ( ftbbe * ftbbe ) - 1);
			//boolean fdkpMass2Cut = (0.2 < FTBmass && FTBmass < 0.40); // cut for K- candidates in FD 
			//boolean cdkpMass2Cut = (0.15 < FTBmass && FTBmass < 0.35); //  cut for K- candidates in CD 
			//boolean kpMass2Cut = (0.180625 < FTBmass && FTBmass < 0.36);// 0.425 < kpM < 0.6
			boolean fdChi2pidCut = (Math.abs(chi2pid) < 5.0);
			boolean cdChi2pidCut = (Math.abs(chi2pid) < 5.0);
			if (pid == 211 ) {
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
				//pip_ftb_beta = ftbbe;
				//pip_status = status;
				
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
								Vpip = new Particle(pid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);
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
								Vpipc = new Particle(pid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);	
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
			if (pid == -211 ) {
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
				//pim_ftb_beta = ftbbe;
				
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
								Vpim = new Particle(pid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);
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
								Vpimc = new Particle(pid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);		
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
			if (pid == 321 ) {
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
				//kp_ftb_beta = ftbbe;
				
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
								Vkp = new Particle(pid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);
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
								Vkpc = new Particle(pid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);
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
			
			if (pid == -321) {
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
				//km_ftb_beta = ftbbe;
				
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
								Vkm = new Particle(pid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);
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
								Vkmc = new Particle(pid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);		
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
			
			  
			 
			
			if (pid == 2212 ) {
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
				//prot_ftb_beta = ftbbe;
				
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
								Vprot = new Particle(pid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);
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
								Vprotc = new Particle(pid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);		
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
					hi_fd_pos_mass.fill(mass);
					H_FD_pos_beta_mom.fill(mom, beta);
					H_FD_pos_mass_mom.fill(mom, mass);
					H_FD_pos_mass_the.fill(the, mass);
					H_FD_pos_mass_phi.fill(phi, mass);		
				}
				if (inCD && cdChi2pidCut) {
					hi_cd_pos_mass.fill(mass);
					H_CD_pos_beta_mom.fill(mom, beta);
					H_CD_pos_mass_mom.fill(mom, mass);
					H_CD_pos_mass_the.fill(the, mass);
					H_CD_pos_mass_phi.fill(phi, mass);				
				}
					
			}
			if ( q < 0 ) {
				if (inDC && fdChi2pidCut) {
					hi_fd_neg_mass.fill(mass);
					H_FD_neg_beta_mom.fill(mom, beta);
					H_FD_neg_mass_mom.fill(mom, mass);
					H_FD_neg_mass_the.fill(the, mass);
					H_FD_neg_mass_phi.fill(phi, mass);		
				}
				if (inCD && cdChi2pidCut) {
					hi_cd_neg_mass.fill(mass);
					H_CD_neg_beta_mom.fill(mom, beta);
					H_CD_neg_mass_mom.fill(mom, mass);
					H_CD_neg_mass_the.fill(the, mass);
					H_CD_neg_mass_phi.fill(phi, mass);				
				}

			}
			if (q == 0 ) {
				if (inDC && fdChi2pidCut) {
					H_FD_neutral_beta_mom.fill(mom, beta);
					H_FD_neutral_mass_mom.fill(mom, mass);
					H_FD_neutral_mass_the.fill(the, mass);
					H_FD_neg_mass_phi.fill(phi, mass);		
				}
				if (inCD && cdChi2pidCut) {
					H_CD_neutral_beta_mom.fill(mom, beta);
					H_CD_neutral_mass_mom.fill(mom, mass);
					H_CD_neutral_mass_the.fill(the, mass);
					H_CD_neutral_mass_phi.fill(phi, mass);				
				}
			}
			
		} // FOR LOOp
		
//		if(npip >= 2 && nfdpip == 2) {
//		System.out.println("found total pip in this event :: " + npip +" with " + nfdpip + " pips tracks in FD " + nfdwithcutpip +" with cut & " + ncdpip + " pip tracks in CD " + ncdwithcutpip + " with cut" );
//		}

	
	} //MAKEOTHER
	
	/*public void makeOthers(DataBank recbank,  DataBank recFTbank, DataBank recSCBank) {
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
	*/

	public boolean select_ekpprot(){
		boolean res = false;
		if(found_eFD && kps.size() == 1 && prots.size() == 1){
			LorentzVector Vmissekpprot = new LorentzVector(0, 0, 0, 0);
			Vmissekpprot.add(VT);Vmissekpprot.add(VB);Vmissekpprot.sub(Ve);Vmissekpprot.sub(kps.get(0).vector());Vmissekpprot.sub(prots.get(0).vector());
			ekpprot_mm2_ekpprot = (float)Vmissekpprot.mass2();

			LorentzVector Vmissekp = new LorentzVector(0, 0, 0, 0);
			Vmissekp.add(VT);Vmissekp.add(VB);Vmissekp.sub(Ve);Vmissekp.sub(kps.get(0).vector());
			ekpprot_mm_ekp = (float)Vmissekp.mass();

			res = true;
		}

		return res;
	}
	

/*
	// function to correct FT electron four vector according to common vertics 
	public Particle electron_corn(Particle el_lv, Vector3D el_FtCal_hitPosition, double vz){

		Vector3D vertex = new Vector3D(el_lv.vx(), el_lv.vy(), vz);
		Vector3D line = new Vector3D(el_FtCal_hitPosition);
		line.sub(vertex);
		//double theta_corn_1 = (el_lv.theta() * vz)/(el_FtCal_hitPosition.z() - vz);
		double thetaCorr = Math.exp(1.797 - 4.485*el_lv.e()) + Math.exp(-0.8671 - 1.078*el_lv.e());
		thetaCorr = Math.toRadians(thetaCorr);
		double phiCorr = Math.exp(4.918 - 3.828*el_lv.e()) + Math.exp(3.841 - 1.256*el_lv.e()) + Math.exp(2.874 - 0.2195*el_lv.e());
		double field = -1; //-1 for negative inbending and +1 for negative outbending. 
		phiCorr = Math.toRadians(phiCorr * field);
		double new_px = (double) el_lv.p()*Math.sin(line.theta() + thetaCorr )*Math.cos(line.phi() - phiCorr);
		double new_py = (double) el_lv.p()*Math.sin(line.theta() + thetaCorr )*Math.sin(line.phi() - phiCorr);
		double new_pz = (double) el_lv.p()*Math.cos(line.theta() + thetaCorr );
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

		//double theta_corn_1 = (el_lv.theta() * vz)/(el_FtCal_hitPosition.z() - vz);
		//double new_px = (double) el_lv.p()*Math.sin(el_lv.theta() + theta_corn_1 )*Math.cos(el_lv.phi());
		//double new_py = (double) el_lv.p()*Math.sin(el_lv.theta() + theta_corn_1)*Math.sin(el_lv.phi());
		//double new_pz = (double) el_lv.p()*Math.cos(el_lv.theta() + theta_corn_1);



		Particle new_Ve = new Particle(11, 0, 0, 0, 0, 0, 0);
		new_Ve.setVector(11, new_px, new_py, new_pz, el_lv.vx(), el_lv.vy(), vz);
		return new_Ve; 
	}

*/

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
		found_eFD = false;
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

		//if (event.hasBank("RECFT::Event"))
		//	fillRecBank(event.getBank("RECFT::Event"));
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
		
		
	} //processEvent
	
	public void analyze() {
		fitxi(dg_rec_y.getH1F("hi_ekpkpkm_MM_ekpkp_nocorr"), dg_rec_y.getF1D("fn_xi_no_corr"));
    	fitxi(dg_rec_y.getH1F("hi_ekpkpkm_MM_ekpkp"), dg_rec_y.getF1D("f1_xi"));
    	fitxi(dg_rec_y.getH1F("hi_mc_ekpkpkm_mm_ekpkp"), dg_rec_y.getF1D("f1_mc_xi"));
    	fitxi(dg_rec_y.getH1F("hi_ekpkpkm_IM_kmlambda"), dg_rec_y.getF1D("fn_im_kmlambda"));
    	fitxi(dg_rec_y.getH1F("hi_ekpkpkm_IM_kmsigma"), dg_rec_y.getF1D("fn_im_kmsigma"));

    	fit_dp(dg_rec_electron.getH1F("hi_rec_e_dp"), dg_rec_electron.getF1D("fn_rec_e_dp"));
    	fit_del(dg_rec_electron.getH1F("hi_rec_e_dtheta"), dg_rec_electron.getF1D("fn_rec_e_dtheta"));
    	fit_del(dg_rec_electron.getH1F("hi_rec_e_dphi"), dg_rec_electron.getF1D("fn_rec_e_dphi"));
    	fit_del(dg_rec_electron.getH1F("hi_rec_e_dvz"), dg_rec_electron.getF1D("fn_rec_e_dvz"));

    	fit_dp(dg_rec_kp.getH1F("hi_rec_kp_dp"), dg_rec_kp.getF1D("fn_rec_kp_dp"));
    	fit_del(dg_rec_kp.getH1F("hi_rec_kp_dtheta"), dg_rec_kp.getF1D("fn_rec_kp_dtheta"));
    	fit_del(dg_rec_kp.getH1F("hi_rec_kp_dphi"), dg_rec_kp.getF1D("fn_rec_kp_dphi"));
    	fit_del(dg_rec_kp.getH1F("hi_rec_kp_dvx"), dg_rec_kp.getF1D("fn_rec_kp_dvx"));
    	fit_del(dg_rec_kp.getH1F("hi_rec_kp_dvy"), dg_rec_kp.getF1D("fn_rec_kp_dvy"));
    	fit_del(dg_rec_kp.getH1F("hi_rec_kp_dvz"), dg_rec_kp.getF1D("fn_rec_kp_dvz"));

    	fit_dp(dg_rec_prot.getH1F("hi_rec_prot_dp"), dg_rec_prot.getF1D("fn_rec_prot_dp"));
    	fit_del(dg_rec_prot.getH1F("hi_rec_prot_dtheta"), dg_rec_prot.getF1D("fn_rec_prot_dtheta"));
    	fit_del(dg_rec_prot.getH1F("hi_rec_prot_dphi"), dg_rec_prot.getF1D("fn_rec_prot_dphi"));
    	fit_del(dg_rec_prot.getH1F("hi_rec_prot_dvx"), dg_rec_prot.getF1D("fn_rec_prot_dvx"));
    	fit_del(dg_rec_prot.getH1F("hi_rec_prot_dvy"), dg_rec_prot.getF1D("fn_rec_prot_dvy"));
    	fit_del(dg_rec_prot.getH1F("hi_rec_prot_dvz"), dg_rec_prot.getF1D("fn_rec_prot_dvz"));

    	fit_dp(dg_rec_km.getH1F("hi_rec_km_dp"), dg_rec_km.getF1D("fn_rec_km_dp"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dtheta"), dg_rec_km.getF1D("fn_rec_km_dtheta"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dphi"), dg_rec_km.getF1D("fn_rec_km_dphi"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dvx"), dg_rec_km.getF1D("fn_rec_km_dvx"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dvy"), dg_rec_km.getF1D("fn_rec_km_dvy"));
    	fit_del(dg_rec_km.getH1F("hi_rec_km_dvz"), dg_rec_km.getF1D("fn_rec_km_dvz"));
    	

    	//fitlambda(dg_rec_y.getH1F("hi_ekpkpkm_MM_ekpkpkm"), dg_rec_y.getF1D("f1_lambda"));
    	//fitlambda(dg_rec_y.getH1F("hi_ekpkpkm_MM_ekpkpkm"), dg_rec_y.getF1D("f1_sigma"));
    	//fitvz(dg_rec_kp.getH1F("hi_rec_kp_dvz"), dg_rec_kp.getF1D("f_gaus"));
    	//fitvz(dg_rec_prot.getH1F("hi_rec_prot_dvz"), dg_rec_prot.getF1D("f_gaus"));
    	//fitvz(dg_rec_km.getH1F("hi_rec_km_dvz"), dg_rec_km.getF1D("f_gaus"));
    //	fitxi(dg_mm2.getH1F("H_"))


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
		
		if(found_eFD){		
//			H_FT_e_t_f.fill(e_phi, e_the);
//			H_FT_e_p_f.fill(e_phi, e_mom);
//			H_FT_e_p_the.fill(e_the, e_mom);
			H_FT_W_Q2.fill(e_W, e_Q2);
			H_FT_W.fill(e_W);
			H_FT_Q2.fill(e_Q2);	
			H_virphoton.fill(e_virphoton);	
		}


		if(select_ekpprot()){

			H_FT_e_t_f.fill(Math.toDegrees(Ve.phi()), Math.toDegrees(Ve.theta()));
			H_FT_e_p_f.fill(Math.toDegrees(Ve.phi()), Ve.p());
			H_FT_e_p_the.fill(Math.toDegrees(Ve.theta()), Ve.p());

			hi_rec_kp_p_the.fill(Math.toDegrees(kps.get(0).theta()), kps.get(0).p());
			hi_rec_kp_p_phi.fill(Math.toDegrees(kps.get(0).phi()), kps.get(0).p());
			hi_rec_kp_the_phi.fill(Math.toDegrees(kps.get(0).phi()), Math.toDegrees(kps.get(0).theta()));
			hi_rec_kp_theta_vz.fill(kps.get(0).vz(), Math.toDegrees(kps.get(0).theta()));
			hi_rec_kp_vz.fill(kps.get(0).vz());

			hi_rec_prot_p_the.fill(Math.toDegrees(prots.get(0).theta()), prots.get(0).p());
			hi_rec_prot_p_phi.fill(Math.toDegrees(prots.get(0).phi()), prots.get(0).p());
			hi_rec_prot_the_phi.fill(Math.toDegrees(prots.get(0).phi()), Math.toDegrees(prots.get(0).theta()));
			hi_rec_prot_theta_vz.fill(prots.get(0).vz(), Math.toDegrees(prots.get(0).theta()));
			hi_rec_prot_vz.fill(prots.get(0).vz());

			hi_ekpprot_mm2_ekpprot.fill(ekpprot_mm2_ekpprot);
			hi_ekpprot_mm_ekp.fill(ekpprot_mm_ekp);
			hi_ekpprot_mm2_ekpprot_mm_ekp.fill(ekpprot_mm_ekp, ekpprot_mm2_ekpprot);


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
		
	public void plot(){

		myCanvas = new EmbeddedCanvasTabbed("Kinematics","REC-electron","REC-kp","REC-prot","REC-Y","FD-TOF", "CD-TOF", "VTime", "particle-Vz","TOF-t","TOF-path", "Counter");
		

		myCanvas.getCanvas("Kinematics").divide(2, 2);
		myCanvas.getCanvas("Kinematics").setGridX(false);
		myCanvas.getCanvas("Kinematics").setGridY(false);
		myCanvas.getCanvas("Kinematics").setAxisFontSize(18);
		myCanvas.getCanvas("Kinematics").setAxisTitleSize(24);
		myCanvas.getCanvas("Kinematics").draw(dg_kinematics);
		myCanvas.getCanvas("Kinematics").getPad(2).getAxisZ().setLog(true);
		//reconstructed electron

		myCanvas.getCanvas("REC-electron").divide(3, 2);
		myCanvas.getCanvas("REC-electron").setGridX(false);
		myCanvas.getCanvas("REC-electron").setGridY(false);
		myCanvas.getCanvas("REC-electron").setAxisFontSize(18);
		myCanvas.getCanvas("REC-electron").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-electron").draw(dg_rec_electron);
		myCanvas.getCanvas("REC-electron").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(2).getAxisZ().setLog(true);

		// reconstructed kp
		myCanvas.getCanvas("REC-kp").divide(3, 2);
		myCanvas.getCanvas("REC-kp").setGridX(false);
		myCanvas.getCanvas("REC-kp").setGridY(false);
		myCanvas.getCanvas("REC-kp").setAxisFontSize(18);
		myCanvas.getCanvas("REC-kp").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-kp").draw(dg_rec_kp);
		myCanvas.getCanvas("REC-kp").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-kp").getPad(2).getAxisZ().setLog(true);
		

		// reconstructed proton
		myCanvas.getCanvas("REC-prot").divide(3, 2);
		myCanvas.getCanvas("REC-prot").setGridX(false);
		myCanvas.getCanvas("REC-prot").setGridY(false);
		myCanvas.getCanvas("REC-prot").setAxisFontSize(18);
		myCanvas.getCanvas("REC-prot").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-prot").draw(dg_rec_prot);
		myCanvas.getCanvas("REC-prot").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-prot").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-prot").getPad(2).getAxisZ().setLog(true);


		myCanvas.getCanvas("REC-Y").divide(2,2);
		myCanvas.getCanvas("REC-Y").setGridX(false);
		myCanvas.getCanvas("REC-Y").setGridY(false);
		myCanvas.getCanvas("REC-Y").setAxisFontSize(18);
		myCanvas.getCanvas("REC-Y").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-Y").draw(dg_rec_y);
		myCanvas.getCanvas("REC-Y").getPad(2).getAxisZ().setLog(true);
		//myCanvas.getCanvas("REC-Xi").getPad(1).draw(L_sigma);

		/*myCanvas.getCanvas("MM2-ekpprot").divide(3,3);
		myCanvas.getCanvas("MM2-ekpprot").setGridX(false);
		myCanvas.getCanvas("MM2-ekpprot").setGridY(false);
		myCanvas.getCanvas("MM2-ekpprot").setAxisFontSize(18);
		myCanvas.getCanvas("MM2-ekpprot").setAxisTitleSize(24);
		myCanvas.getCanvas("MM2-ekpprot").draw(dg_mm2);
		myCanvas.getCanvas("MM2-ekpprot").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM2-ekpprot").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM2-ekpprot").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM2-ekpprot").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM2-ekpprot").getPad(4).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM2-ekpprot").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM2-ekpprot").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM2-ekpprot").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM2-ekpprot").getPad(8).getAxisZ().setLog(true);
*/
		// TOF
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
		efdky_data_v1 ana = new efdky_data_v1();
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
		//ana.analyze();
		ana.plot();
		ana.showplots();
		//ana.save();
		System.out.println("Good Bye !!!!!");
	}
		
	
	
}
