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
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
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
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
//import org.jlab.io.hipo.*;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import javax.swing.JFrame;

//import from hayward_coatjava_extensions 
import org.jlab.clas.physics.GenericKinematicFitter;
import extended_kinematic_fitters.*;

/**
 * @author akhanal
 *
 */

public class pionsaskaons_ekpkpreqkmavgvz_data_v1 {
	//public int NFTElec;

	//public analysis_fitter fiducial_cut;
	
	public static float  Eb, Mp;
	public static LorentzVector VB;
	public static LorentzVector VT;

	//Timothy's analysis fitter
	public static GenericKinematicFitter research_fitter;
	public static File outFile, me_outFile, me_twoKaon_outFile;
	public static File mcoutFile;
	public static String runType;

	public static ArrayList<LorentzVector> lv_my_ftels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
	public static ArrayList<LorentzVector> lv_my_xi_ftels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms;

	//public float EB, Eb, Mp;
	public float STT, RFT, FTSTT, vt;
	public static int runnum, evnum;
	public static double smearFactor;

	public Particle Vprot, Vpip, Vpim, Vkp, Vkm;
	public Particle Vprotc, Vpipc, Vpimc, Vkpc, Vkmc;
	public Particle Velectron, Vftel_corrected;
	public LorentzVector Ve_consCorn;
	public LorentzVector Ve, VGS, VhadronSystm, Vpim_correct, Vkp_correct;

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
	
	
	public H1F hi_FT_W, hi_FT_Q2, hi_virphoton;
	public H1F  hi_ekpkp_MM_req_pim, hi_efkpckp_MM_req_pim, hi_efkpckp_MM_req_cpim, hi_efkpckp_MM_req_km, hi_efkpckp_MM_req_ckm;
	// particles vz, pathi_FTOF1b, time_FTOF1b
	
	public H1F hi_kp_vz, hi_km_vz, hi_prot_vz, hi_pip_vz, hi_pim_vz;
	public H1F hi_pipc_vz, hi_pimc_vz, hi_kpc_vz, hi_kmc_vz, hi_protc_vz; 
	public H1F hi_prot_FTOF1b_path, hi_pip_FTOF1b_path,  hi_pim_FTOF1b_path, hi_kp_FTOF1b_path, hi_km_FTOF1b_path;
	public H1F hi_prot_FTOF1b_t, hi_pip_FTOF1b_t,  hi_pim_FTOF1b_t, hi_kp_FTOF1b_t, hi_km_FTOF1b_t;
	public H1F hi_pip_CTOF_path, hi_pim_CTOF_path, hi_prot_CTOF_path,  hi_kp_CTOF_path, hi_km_CTOF_path;
	public H1F hi_pim_CTOF_t, hi_pip_CTOF_t, hi_prot_CTOF_t, hi_kp_CTOF_t, hi_km_CTOF_t;
	

	public H2F hi_pip_vt_p, hi_prot_vt_p, hi_pim_vt_p, hi_kp_vt_p, hi_km_vt_p;
	public H2F hi_pipc_vt_p, hi_pimc_vt_p, hi_kpc_vt_p, hi_kmc_vt_p, hi_protc_vt_p;

	public H1F hi_cd_prot_p, hi_cd_pip_p;
	public H1F hi_cd_prot_theta, hi_cd_pip_theta;
	public H2F hi_cd_prot_p_theta, hi_cd_pip_p_theta;
	public H2F hi_W_cd_pro_the;

	
	public H1F hi_fd_pos_mass, hi_fd_neg_mass, hi_cd_pos_mass, hi_cd_neg_mass;
	public H2F hi_FD_pos_beta_mom, hi_FD_neg_beta_mom, hi_FD_neutral_beta_mom;
	public H2F hi_FD_pos_mass_mom, hi_FD_neg_mass_mom, hi_FD_neutral_mass_mom;
	public H2F hi_FD_pos_mass_the, hi_FD_neg_mass_the, hi_FD_neutral_mass_the;
	public H2F hi_FD_pos_mass_phi, hi_FD_neg_mass_phi, hi_FD_neutral_mass_phi;
	public H2F hi_CD_pos_beta_mom, hi_CD_neg_beta_mom, hi_CD_neutral_beta_mom;
	public H2F hi_CD_pos_mass_mom, hi_CD_neg_mass_mom, hi_CD_neutral_mass_mom;
	public H2F hi_CD_pos_mass_the, hi_CD_neg_mass_the, hi_CD_neutral_mass_the;
	public H2F hi_CD_pos_mass_phi, hi_CD_neg_mass_phi, hi_CD_neutral_mass_phi;
	/*
	public List<Particle> fpips, fpims, fkps, fkms, fprots; // = new ArrayList<Particle>();
	public List<Particle> cpips, cpims, ckps, ckms, cprots;
	public List<Particle> ftels, pips, pims, kps, kms, prots;
	public List<Particle> mcels, mckps, mckms;
	*/

	public ArrayList<Particle> fpips, fpims, fkps, fkms, fprots; // = new ArrayList<Particle>();
	public ArrayList<Particle> cpips, cpims, ckps, ckms, cprots;
	public ArrayList<Particle> ftels, pips, pims, kps, kms, prots;
	public ArrayList<Particle> mcels, mckps, mckms;

    public Particle mckp1;
    public Particle mckp2;
	// for e kp kp km detected
	public Particle reckp1;
    public Particle reckp2;
	public float rec_kp1_p, rec_kp1_the, rec_kp1_phi, rec_kp1_vz, rec_kp2_p, rec_kp2_the, rec_kp2_phi, rec_kp2_vz;
	public float rec_km_p, rec_km_the, rec_km_phi, rec_km_vz;
	public float efkpfkpfkm_MM_efkpfkp,  efkpckpfkm_MM_efkpckp,  eckpckpfkm_MM_eckpckp;  
	public float efkpfkpckm_MM_efkpfkp,  efkpckpckm_MM_efkpckp,  eckpckpckm_MM_eckpckp;
	public float efkpfkpfkm_MM_efkpfkpfkm,  efkpckpfkm_MM_efkpckpfkm, efkpckpfkm_MM_efkpckpfkm_nocorn,  eckpckpfkm_MM_eckpckpfkm;  
	public float efkpfkpckm_MM_efkpfkpckm,  efkpckpckm_MM_efkpckpckm,  eckpckpckm_MM_eckpckpckm;
	public float efkpfkpfkm_IM_kmlambda, efkpckpfkm_IM_kmlambda, eckpckpfkm_IM_kmlambda;
	public float efkpfkpckm_IM_kmlambda, efkpckpckm_IM_kmlambda, eckpckpckm_IM_kmlambda;
	public float ekpkpkm_MM_ekpkp, ekpkpkm_MM_ekpkpkm, ekpkpkm_IM_kmlambda, ekpkpkm_IM_kmsigma;
	public float ekpkpkm_MM_ekpkp_nocorr, ekpkpkm_MM_ekpkpkm_nocorr;
	public float ekpkpkmprot_MM2;

	/// various correction functions for particle in CD 
	public F1D fn_deltaphi_phi1, fn_deltaphi_phi2, fn_deltaphi_phi3, fn_deltaphi_theta, fn_deltaphi_p;
    public F1D fn_deltatheta_phi1, fn_deltatheta_phi2, fn_deltatheta_phi3, fn_deltatheta_theta, fn_deltatheta_p;
    public F1D fn_deltap_phi1, fn_deltap_phi2, fn_deltap_phi3, fn_deltap_theta, fn_deltap_p;

    public F1D fn2_deltaphi_phi1, fn2_deltaphi_phi2, fn2_deltaphi_phi3;
    public F1D fn2_deltatheta_theta, fn2_deltatheta_p;
	
	public F1D fn_rec_e_dp, fn_rec_e_dtheta, fn_rec_e_dphi, fn_rec_e_dvx, fn_rec_e_dvy, fn_rec_e_dvz;
	public F1D fn_rec_kp1_dp, fn_rec_kp1_dtheta, fn_rec_kp1_dphi, fn_rec_kp1_dvx, fn_rec_kp1_dvy, fn_rec_kp1_dvz;
	public F1D fn_rec_kp2_dp, fn_rec_kp2_dtheta, fn_rec_kp2_dphi, fn_rec_kp2_dvx, fn_rec_kp2_dvy, fn_rec_kp2_dvz;
	public F1D fn_rec_km_dp, fn_rec_km_dtheta, fn_rec_km_dphi, fn_rec_km_dvx, fn_rec_km_dvy, fn_rec_km_dvz;

	public H2F hi_FT_e_beta_mom;
	public H2F hi_FT_e_t_f, hi_FT_e_p_f, hi_FT_e_p_the;
	public H1F hi_rec_e_dp, hi_rec_e_dtheta, hi_rec_e_dphi;
	public H2F hi_rec_e_dp_p, hi_rec_e_dp_theta, hi_rec_e_dp_phi, hi_rec_e_dp_vz, hi_rec_e_dtheta_p, hi_rec_e_dtheta_theta, hi_rec_e_dtheta_phi, hi_rec_e_dtheta_vz, hi_rec_e_dphi_p, hi_rec_e_dphi_theta, hi_rec_e_dphi_phi, hi_rec_e_dphi_vz, hi_rec_e_dvz_p, hi_rec_e_dvz_theta, hi_rec_e_dvz_phi, hi_rec_e_dvz_vz;
	public H2F hi_FT_W_Q2, hi_FT_e_xB_Q2;

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
	public H1F hi_ekpkpkm_MM_ekpkp_nocorr, hi_ekpkpkm_MM_ekpkpkm_nocorr
	public F1D f1_xi, fn_xi_no_corr, f1_mc_xi, f1_lambda, f1_sigma, f_gaus, fn_im_kmlambda, fn_im_kmsigma;
	public H1F hi_mc_ekpkpkm_mm_ekpkp, hi_mc_ekpkpkm_mm_ekpkpkm;

	//mass spectrum ekpkpkm detected
	public H1F hi_efkpfkpfkm_MM_efkpfkp, hi_efkpckpfkm_MM_efkpckp, hi_eckpckpfkm_MM_eckpckp;
	public H1F hi_efkpfkpckm_MM_efkpfkp, hi_efkpckpckm_MM_efkpckp, hi_eckpckpckm_MM_eckpckp;

	public H1F hi_fkpckp_deltap, hi_fkpckp_deltap_withdpdtcut, hi_fkpckp_deltatheta, hi_fkpckp_deltatheta_withdpdtcut, hi_fkpckp_deltaphi, hi_fkpckp_deltaphi_withdpdtcut;
	public H2F hi_fkpckp_deltap_p, hi_fkpckp_deltap_theta, hi_fkpckp_deltap_phi;
	public H2F hi_fkpckp_deltatheta_p, hi_fkpckp_deltatheta_theta, hi_fkpckp_deltatheta_phi;
	public H2F hi_fkpckp_deltaphi_p, hi_fkpckp_deltaphi_theta, hi_fkpckp_deltaphi_phi;

	public H1F hi_efkpfkpfkm_MM_efkpfkpfkm, hi_efkpckpfkm_MM_efkpckpfkm, hi_efkpckpfkm_MM_efkpckpfkm_nocorn, hi_eckpckpfkm_MM_eckpckpfkm;
	public H1F hi_efkpfkpckm_MM_efkpfkpckm, hi_efkpckpckm_MM_efkpckpckm, hi_eckpckpckm_MM_eckpckpckm;

	//mass spectra with cut(lambda event in MM(ekpkpkm))/ lambda constraint
	public H1F hi_efkpfkpfkm_MM_efkpfkp_lam_evnt, hi_efkpckpfkm_MM_efkpckp_lam_evnt, hi_eckpckpfkm_MM_eckpckp_lam_evnt;
	public H1F hi_efkpfkpckm_MM_efkpfkp_lam_evnt, hi_efkpckpckm_MM_efkpckp_lam_evnt, hi_eckpckpckm_MM_eckpckp_lam_evnt;
	public H1F hi_efkpfkpfkm_IM_kmlambda, hi_efkpckpfkm_IM_kmlambda, hi_eckpckpfkm_IM_kmlambda;
	public H1F hi_efkpfkpckm_IM_kmlambda, hi_efkpckpckm_IM_kmlambda, hi_eckpckpckm_IM_kmlambda;


	public boolean efkpfkpfkm_found_lambda, efkpckpfkm_found_lambda, eckpckpfkm_found_lambda;
	public boolean efkpfkpckm_found_lambda, efkpckpckm_found_lambda, eckpckpckm_found_lambda;


	// scatter plots
	public H2F hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm, hi_efkpckpfkm_MM_efkpckp_MM_efkpckpfkm, hi_eckpckpfkm_MM_eckpckp_MM_eckpckpfkm;
	public H2F hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm, hi_efkpckpckm_MM_efkpckp_MM_efkpckpckm, hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm;


	public H1F hi_ekpkpkmprot_MM2;

	public H1F hi_pip_counter, hi_pim_counter, hi_kp_counter, hi_km_counter, hi_prot_counter, hi_fpip_counter, hi_fpim_counter, hi_fkp_counter, hi_fkm_counter, hi_fprot_counter, hi_cpip_counter, hi_cpim_counter, hi_ckp_counter, hi_ckm_counter, hi_cprot_counter;

	//public EmbeddedCanvas rec_electron;
	public EmbeddedCanvasTabbed myCanvas;
	//public EmbeddedCanvas myCanvas;
	public DataGroup dg_rec_electron, dg_rec_kp1, dg_rec_kp2, dg_rec_km, dg_rec_xi, dg_rec_p, dg_rec_pim, dg_vtime, dg_fdtof, dg_cdtof, dg_cdPart, dg_vz, dg_tof_t, dg_tof_path, dg_counter;
	public DataGroup dg_rec_e_resolution, dg_rec_kp1_resolution, dg_rec_kp2_resolution, dg_rec_km_resolution;
	public DataGroup dg_mm, dg_mm_ekpkpkm, dg_mm_ekpkpkm1, dg_mm_scatter, dg_mm_ekpkpkm_lam_evnt, dg_m_kmlambda, dg_fkpckp_pthphi;
	public DataGroup dg_req;

	
	
	public pionsaskaons_ekpkpreqkmavgvz_data_v1() {

		final int RED = 2;
		final int BLUE = 9;
		final int LIGHTGREEN = 3;
		final int LIGHTBROWN = 45;
		final int PINK = 46;
	//	NFTElec = 0;
	//	Eb = 10.575f;
	//	Eb = 10.604f; //RGA fall2018 beam energy
	//	Eb = 7.54626f;
	//	Eb = 10.1998; //RGA spring2019 beam energy
		Mp = (float) PDGDatabase.getParticleById(2212).mass();
	//	Mp = 0.93827f;
		
	//	VB = new LorentzVector(0, 0, Eb, Eb);
		VT = new LorentzVector(0, 0, 0, Mp);
	
		
		// theoretical 1D functions for proton, kaon and pion
		
		//F_prot_beta_mom = new F1D("F_prot_beta_mom", "x/sqrt(0.93827*0.93827+x*x)", 0.3, 4.0);
		F_prot_beta_mom = new F1D("F-prot-beta-mom", "x/sqrt(0.93827*0.93827+x*x)", 0.1, 5.0);
		F_prot_beta_mom.setLineWidth(2);
		F_prot_beta_mom.setLineColor(BLUE);
		//F_kp_beta_mom = new F1D("F_kp_beta_mom", "x/sqrt(0.49367*0.49367+x*x)", 0.3, 4.0);
		F_kp_beta_mom = new F1D("F-kp-beta-mom", "x/sqrt(0.49367*0.49367+x*x)", 0.1, 5.0);
		F_kp_beta_mom.setLineColor(BLUE);
		F_kp_beta_mom.setLineWidth(2);
		//F_pip_beta_mom = new F1D("F_pip_beta_mom", "x/sqrt(0.13957*0.13957+x*x)", 0.3, 4.0);
		F_pip_beta_mom = new F1D("F-pip-beta-mom", "x/sqrt(0.13957*0.13957+x*x)", 0.1, 5.0);
		F_pip_beta_mom.setLineColor(BLUE);
		F_pip_beta_mom.setLineWidth(2);


		// CD particle correction functions derived using central pip in fefpcpipfpim events considering cpip is missing

		// phi correction functions for three regions of CVT from deltaPhi Vs phi plot using cpip in eppippim event
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

        //momentum correction functions for three regions of CVT from deltatheta Vs phi plot
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
		hi_FT_e_t_f = new H2F("hi_FT_e_t_f", "hi_FT_e_t_f", 100, -180, 180, 100, 2, 5);
		hi_FT_e_t_f.setTitle("electron #theta vs #phi");
		hi_FT_e_t_f.setTitleX("#phi (^o)");
		hi_FT_e_t_f.setTitleY("#theta (^o)");
		
		hi_FT_e_p_the = new H2F("hi_FT_e_p_the", "hi_FT_e_p_the", 100, 2, 5, 100, 0, 5);
		hi_FT_e_p_the.setTitle("electron p vs #theta (^o)");
		hi_FT_e_p_the.setTitleX("#theta (^o)");
		hi_FT_e_p_the.setTitleY("p (GeV)");

		hi_FT_e_p_f = new H2F("hi_FT_e_p_f", "hi_FT_e_p_f", 100, -180, 180, 100, 0, 5);
		hi_FT_e_p_f.setTitle("electron p vs #phi");
		hi_FT_e_p_f.setTitleX("#phi (^o)");
		hi_FT_e_p_f.setTitleY("p (GeV)");

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
		
		hi_FT_W_Q2 = new H2F("hi_FT_W_Q2", "hi_FT_W_Q2", 100, 3.45, 4.55, 100, 0.005, 0.25);
		hi_FT_W_Q2.setTitle("FT Q^2 vs W");
		hi_FT_W_Q2.setTitleX("W ( GeV)");
		hi_FT_W_Q2.setTitleY("Q^2 (GeV^2)");
		
		
		hi_FT_W = new H1F("hi_FT_W", "hi_FT_W", 100, 3, 5);
		hi_FT_W.setTitle("electron W");
		hi_FT_W.setTitleX("W (GeV)");
		hi_FT_W.setTitleY("count");
		hi_FT_W.setFillColor(LIGHTGREEN);


		hi_FT_Q2 = new H1F("hi_FT_Q2", "hi_FT_Q2", 100, 0.00001, 0.5);
		hi_FT_Q2.setFillColor(LIGHTGREEN);
		hi_FT_Q2.setTitleX("Q^2 (GeV^2)");

		hi_virphoton = new H1F("hi_virphoton", "hi_virphoton", 100, 0, 12);
		hi_virphoton.setFillColor(LIGHTGREEN);
		hi_virphoton.setTitleX("E_#gamma (GeV)");

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

		//ekpkpkm required
		hi_efkpfkpfkm_MM_efkpfkp = new H1F("hi-efkpfkpfkm-MM-efkpfkp", "hi-efkpfkpfkm-MM-efkpfkp", 50, 1.6, 3.1); //50, 1.6, 2.1 //100, 1.6, 3.1//75, 1.6, 3.1
		hi_efkpfkpfkm_MM_efkpfkp.setFillColor(LIGHTGREEN);
		hi_efkpckpfkm_MM_efkpckp = new H1F("hi-efkpckpfkm-MM-efkpckp", "hi-efkpckpfkm-MM-efkpckp", 50, 1.6, 3.1);
		hi_efkpckpfkm_MM_efkpckp.setFillColor(LIGHTGREEN);
		hi_eckpckpfkm_MM_eckpckp = new H1F("hi-eckpckpfkm-MM-eckpckp", "hi-eckpckpfkm-MM-eckpckp", 50, 1.6, 3.1);
		hi_eckpckpfkm_MM_eckpckp.setFillColor(LIGHTGREEN);
		hi_efkpfkpckm_MM_efkpfkp = new H1F("hi-efkpfkpckm-MM-efkpfkp", "hi-efkpfkpckm-MM-efkpfkp", 50, 1.6, 3.1);
		hi_efkpfkpckm_MM_efkpfkp.setFillColor(LIGHTGREEN);
		hi_efkpckpckm_MM_efkpckp = new H1F("hi-efkpckpckm-MM-efkpckp", "hi-efkpckpckm-MM-efkpckp", 50, 1.6, 3.1);
		hi_efkpckpckm_MM_efkpckp.setFillColor(LIGHTGREEN);
		hi_eckpckpckm_MM_eckpckp = new H1F("hi-eckpckpckm-MM-eckpckp", "hi-eckpckpckm-MM-eckpckp", 50, 1.6, 3.1);
		hi_eckpckpckm_MM_eckpckp.setFillColor(LIGHTGREEN);


		// with lambda events selection in MM(ekpkpkm)
		hi_efkpfkpfkm_MM_efkpfkp_lam_evnt = new H1F("hi-efkpfkpfkm-MM-efkpfkp-lam-evnt", "hi-efkpfkpfkm-MM-efkpfkp-lam-evnt", 30, 1.6, 2.5);// 30, 1.6, 2.5 //25, 1.6, 2.6 //50, 1.6, 2.1 //100, 1.6, 3.1//75, 1.6, 3.1
		hi_efkpfkpfkm_MM_efkpfkp_lam_evnt.setFillColor(LIGHTGREEN);
		hi_efkpckpfkm_MM_efkpckp_lam_evnt= new H1F("hi-efkpckpfkm-MM-efkpckp-lam-evnt", "hi-efkpckpfkm-MM-efkpckp-lam-evnt", 30, 1.6, 2.5);
		hi_efkpckpfkm_MM_efkpckp_lam_evnt.setFillColor(LIGHTGREEN);
		hi_eckpckpfkm_MM_eckpckp_lam_evnt = new H1F("hi-eckpckpfkm-MM-eckpckp-lam-evnt", "hi-eckpckpfkm-MM-eckpckp-lam-evnt", 30, 1.6, 2.5);
		hi_eckpckpfkm_MM_eckpckp_lam_evnt.setFillColor(LIGHTGREEN);
		hi_efkpfkpckm_MM_efkpfkp_lam_evnt = new H1F("hi-efkpfkpckm-MM-efkpfkp-lam-evnt", "hi-efkpfkpckm-MM-efkpfkp-lam-evnt", 30, 1.6, 2.5);
		hi_efkpfkpckm_MM_efkpfkp_lam_evnt.setFillColor(LIGHTGREEN);
		hi_efkpckpckm_MM_efkpckp_lam_evnt = new H1F("hi-efkpckpckm-MM-efkpckp-lam-evnt", "hi-efkpckpckm-MM-efkpckp-lam-evnt", 30, 1.6, 2.5);
		hi_efkpckpckm_MM_efkpckp_lam_evnt.setFillColor(LIGHTGREEN);
		hi_eckpckpckm_MM_eckpckp_lam_evnt = new H1F("hi-eckpckpckm-MM-eckpckp-lam-evnt", "hi-eckpckpckm-MM-eckpckp-lam-evnt", 30, 1.6, 2.5);
		hi_eckpckpckm_MM_eckpckp_lam_evnt.setFillColor(LIGHTGREEN);

		// IM (kmlambdaconstrained) 
		hi_efkpfkpfkm_IM_kmlambda = new H1F("hi-efkpfkpfkm-IM-kmlambda", "hi-efkpfkpfkm-IM-kmlambda", 30, 1.6, 2.5);//15, 1.6, 2.5//30, 1.6, 2.5 //30, 1.7, 2.6//33, 1.6, 2.6//25, 1.6, 2.6//50, 1.6, 3.1 //50, 1.6, 2.1 //100, 1.6, 3.1//75, 1.6, 3.1
		hi_efkpfkpfkm_IM_kmlambda.setFillColor(LIGHTGREEN);
		hi_efkpckpfkm_IM_kmlambda = new H1F("hi-efkpckpfkm-IM-kmlambda", "hi-efkpckpfkm-IM-kmlambda", 30, 1.6, 2.5);
		hi_efkpckpfkm_IM_kmlambda.setFillColor(LIGHTGREEN);
		hi_eckpckpfkm_IM_kmlambda = new H1F("hi-eckpckpfkm-IM-kmlambda", "hi-eckpckpfkm-IM-kmlambda", 30, 1.6, 2.5);
		hi_eckpckpfkm_IM_kmlambda.setFillColor(LIGHTGREEN);
		hi_efkpfkpckm_IM_kmlambda = new H1F("hi-efkpfkpckm-IM-kmlambda", "hi-efkpfkpckm-IM-kmlambda", 30, 1.6, 2.5);
		hi_efkpfkpckm_IM_kmlambda.setFillColor(LIGHTGREEN);
		hi_efkpckpckm_IM_kmlambda = new H1F("hi-efkpckpckm-IM-kmlambda", "hi-efkpckpckm-IM-kmlambda", 30, 1.6, 2.5);
		hi_efkpckpckm_IM_kmlambda.setFillColor(LIGHTGREEN);
		hi_eckpckpckm_IM_kmlambda = new H1F("hi-eckpckpckm-IM-kmlambda", "hi-eckpckpckm-IM-kmlambda", 30, 1.6, 2.5);
		hi_eckpckpckm_IM_kmlambda.setFillColor(LIGHTGREEN);

		hi_fkpckp_deltap = new H1F("hi-fkpckp-deltap","hi-fkpckp-deltap",100, -1, 1);
		hi_fkpckp_deltap.setTitleX("#Deltap");
		hi_fkpckp_deltap.setFillColor(LIGHTGREEN);
		hi_fkpckp_deltap_withdpdtcut = new H1F("hi-fkpckp-deltap-withdpdtcut","hi-fkpckp-deltap-withdpdtcut",100, -1, 1);
		hi_fkpckp_deltap_withdpdtcut.setFillColor(LIGHTGREEN);
		hi_fkpckp_deltap_withdpdtcut.setLineColor(RED);
		hi_fkpckp_deltatheta = new H1F("hi-fkpckp-deltatheta","hi-fkpckp-deltatheta",100, -100, 40);
		hi_fkpckp_deltatheta.setTitleX("#Delta#theta");
		hi_fkpckp_deltatheta.setFillColor(LIGHTGREEN);
		hi_fkpckp_deltatheta_withdpdtcut = new H1F("hi-fkpckp-deltatheta-withdpdtcut","hi-fkpckp-deltatheta-withdpdtcut", 100, -100, 40);
		hi_fkpckp_deltatheta_withdpdtcut.setFillColor(LIGHTGREEN);
		hi_fkpckp_deltatheta_withdpdtcut.setLineColor(RED);
		hi_fkpckp_deltaphi = new H1F("hi-fkpckp-deltaphi","hi-fkpckp-deltaphi",100, -100, 100);
		hi_fkpckp_deltaphi.setTitleX("#Delta#phi");
		hi_fkpckp_deltaphi.setFillColor(LIGHTGREEN);
		hi_fkpckp_deltaphi_withdpdtcut = new H1F("hi-fkpckp-deltaphi-withdpdtcut","hi-fkpckp-deltaphi-withdpdtcut",100, -100, 100);
		hi_fkpckp_deltaphi_withdpdtcut.setFillColor(LIGHTGREEN);
		hi_fkpckp_deltaphi_withdpdtcut.setLineColor(RED);

		hi_fkpckp_deltap_p = new H2F("hi-fkpckp-deltap-p", "hi-fkpckp-deltap-p", 10, 0.4, 1.8, 30, -0.45, 0.1);
		hi_fkpckp_deltap_p.setTitleX("p");
		hi_fkpckp_deltap_p.setTitleY("#DeltaP");
		hi_fkpckp_deltap_theta = new H2F("hi-fkpckp-deltap-theta", "hi-fkpckp-deltap-theta", 10, 34, 46, 30, -0.45, 0.1);
		hi_fkpckp_deltap_theta.setTitleX("#theta");
		hi_fkpckp_deltap_theta.setTitleY("#DeltaP");
		hi_fkpckp_deltap_phi = new H2F("hi-fkpckp-deltap-phi", "hi-fkpckp-deltap-phi", 15, -180, 180, 30, -0.45, 0.1);
		hi_fkpckp_deltap_phi.setTitleX("#phi");
		hi_fkpckp_deltap_phi.setTitleY("#DeltaP");

		hi_fkpckp_deltatheta_p = new H2F("hi-fkpckp-deltatheta-p", "hi-fkpckp-deltatheta-p", 10, 0.4, 1.8, 30, -8, 8);
		hi_fkpckp_deltatheta_p.setTitleX("p");
		hi_fkpckp_deltatheta_p.setTitleY("#Delta#theta");
		hi_fkpckp_deltatheta_theta = new H2F("hi-fkpckp-deltatheta-theta", "hi-fkpckp-deltatheta-theta", 10, 34, 46, 30, -8, 8);
		hi_fkpckp_deltatheta_theta.setTitleX("#theta");
		hi_fkpckp_deltatheta_theta.setTitleY("#Delta#theta");
		hi_fkpckp_deltatheta_phi = new H2F("hi-fkpckp-deltatheta-phi", "hi-fkpckp-deltatheta-phi", 15, -180, 180, 30, -8, 8);
		hi_fkpckp_deltatheta_phi.setTitleX("#phi");
		hi_fkpckp_deltatheta_phi.setTitleY("#Delta#theta");

		hi_fkpckp_deltaphi_p = new H2F("hi-fkpckp-deltaphi-p", "hi-fkpckp-deltaphi-p", 10, 0.4, 1.8, 30, -25, 10);
		hi_fkpckp_deltaphi_p.setTitleX("p");
		hi_fkpckp_deltaphi_p.setTitleY("#Delta#phi");
		hi_fkpckp_deltaphi_theta = new H2F("hi-fkpckp-deltaphi-theta", "hi-fkpckp-deltaphi-theta", 10, 34, 46, 30, -25, 10);
		hi_fkpckp_deltaphi_theta.setTitleX("#theta");
		hi_fkpckp_deltaphi_theta.setTitleY("#Delta#phi");
		hi_fkpckp_deltaphi_phi = new H2F("hi-fkpckp-deltaphi-phi", "hi-fkpckp-deltaphi-phi", 15, -180, 180, 30, -25, 10);
		hi_fkpckp_deltaphi_phi.setTitleX("#phi");
		hi_fkpckp_deltaphi_phi.setTitleY("#Delta#phi");

		hi_efkpfkpfkm_MM_efkpfkpfkm = new H1F("hi-efkpfkpfkm-MM-efkpfkpfkm", "hi-efkpfkpfkm-MM-efkpfkpfkm", 100, 0.0, 2.4); //25, 0.0, 2.4//25, 0.9, 1.4 //400, 0, 4
		hi_efkpfkpfkm_MM_efkpfkpfkm.setFillColor(LIGHTGREEN);
		hi_efkpckpfkm_MM_efkpckpfkm = new H1F("hi-efkpckpfkm-MM-efkpckpfkm", "hi-efkpckpfkm-MM-efkpckpfkm", 100, 0.0, 2.4); //25, 0.0, 2.4//25, 0.9, 1.4//100, 0.0, 2.4
		hi_efkpckpfkm_MM_efkpckpfkm.setLineColor(RED);
		hi_efkpckpfkm_MM_efkpckpfkm_nocorn = new H1F("hi-efkpckpfkm-MM-efkpckpfkm-nocorn", "hi-efkpckpfkm-MM-efkpckpfkm-nocorn", 100, 0.0, 2.4);
		hi_efkpckpfkm_MM_efkpckpfkm_nocorn.setFillColor(LIGHTGREEN);
		hi_eckpckpfkm_MM_eckpckpfkm = new H1F("hi-eckpckpfkm-MM-eckpckpfkm", "hi-eckpckpfkm-MM-eckpckpfkm", 25, 0.9, 1.4);
		hi_eckpckpfkm_MM_eckpckpfkm.setFillColor(LIGHTGREEN);
		hi_efkpfkpckm_MM_efkpfkpckm = new H1F("hi-efkpfkpckm-MM-efkpfkpckm", "hi-efkpfkpckm-MM-efkpfkpckm", 25, 0.9, 1.4);
		hi_efkpfkpckm_MM_efkpfkpckm.setFillColor(LIGHTGREEN);
		hi_efkpckpckm_MM_efkpckpckm = new H1F("hi-efkpckpckm-MM-efkpckpckm", "hi-efkpckpckm-MM-efkpckpckm", 25, 0.9, 1.4);
		hi_efkpckpckm_MM_efkpckpckm.setFillColor(LIGHTGREEN);
		hi_eckpckpckm_MM_eckpckpckm = new H1F("hi-eckpckpckm-MM-eckpckpckm", "hi-eckpckpckm-MM-eckpckpckm", 25, 0.9, 1.4);
		hi_eckpckpckm_MM_eckpckp.setFillColor(LIGHTGREEN);

		// scatter plots
		hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm = new H2F("hi-efkpfkpfkm-MM-efkpfkp-MM-efkpfkpfkm", "hi-efkpfkpfkm-MM-efkpfkp-MM-efkpfkpfkm", 25, 0.9, 1.4, 25, 1.6, 3.1); //50, 1.6, 2.1 //100, 1.6, 3.1//75, 1.6, 3.1
		//hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm.setFillColor(LIGHTGREEN);
		hi_efkpckpfkm_MM_efkpckp_MM_efkpckpfkm = new H2F("hi-efkpckpfkm-MM-efkpckp-MM-efkpckpfkm", "hi-efkpckpfkm-MM-efkpckp-MM-efkpckpfkm", 25, 0.9, 1.4, 25, 1.6, 3.1);
		//hi_efkpckpfkm_MM_efkpckp_MM_efkpckpfkm.setFillColor(LIGHTGREEN);
		hi_eckpckpfkm_MM_eckpckp_MM_eckpckpfkm = new H2F("hi-eckpckpfkm-MM-eckpckp-MM-eckpckpfkm", "hi-eckpckpfkm-MM-eckpckp-MM-eckpckpfkm", 25, 0.9, 1.4, 25, 1.6, 3.1);
		//hi_eckpckpfkm_MM_eckpckp_MM_eckpckpfkm.setFillColor(LIGHTGREEN);
		hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm = new H2F("hi-efkpfkpckm-MM-efkpfkp-MM-efkpfkpckm", "hi-efkpfkpckm-MM-efkpfkp-MM-efkpfkpckm", 25, 0.9, 1.4, 25, 1.6, 3.1);
		//hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm.setFillColor(LIGHTGREEN);
		hi_efkpckpckm_MM_efkpckp_MM_efkpckpckm = new H2F("hi-efkpckpckm-MM-efkpckp-MM-efkpckpckm", "hi-efkpckpckm-MM-efkpckp-MM-efkpckpckm", 25, 0.9, 1.4, 25, 1.6, 3.1);
		//hi_efkpckpckm_MM_efkpckp_MM_efkpckpckm.setFillColor(LIGHTGREEN);
		hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm = new H2F("hi-eckpckpckm-MM-eckpckp-MM-eckpckpckm", "hi-eckpckpckm-MM-eckpckp-MM-eckpckpckm", 25, 0.9, 1.4, 25, 1.6, 3.1);
		//hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm.setFillColor(LIGHTGREEN);
		
	//ekpkpreqkm	reconstructed
		
		hi_ekpkpkm_MM_ekpkp = new H1F("hi_ekpkpkm_MM_ekpkp", "hi_ekpkpkm_MM_ekpkp", 25, 1.6, 2.1);
		hi_ekpkpkm_MM_ekpkp.setTitle("MM");
		hi_ekpkpkm_MM_ekpkp.setTitleX("MM(eK^+K^+) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkp.setTitleY("Events/[20 MeV/c^2]");
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
		hi_ekpkpkm_MM_ekpkpkm = new H1F("hi_ekpkpkm_MM_ekpkpkm", "hi_ekpkpkm_MM_ekpkpkm", 25, 0.9, 1.4);
		hi_ekpkpkm_MM_ekpkpkm.setTitle("MM");
		hi_ekpkpkm_MM_ekpkpkm.setTitleX("MM(eK^+K^+K^-) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm.setTitleY("Events/[20 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm.setFillColor(LIGHTGREEN);

		// for one with common vertex correction for electron
		hi_ekpkpkm_MM_ekpkpkm_nocorr = new H1F("hi_ekpkpkm_MM_ekpkpkm_nocorr", "hi_ekpkpkm_MM_ekpkpkm_nocorr", 25, 0.9, 1.4);
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitle("MM");
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitleX("MM(eK^+K^+K^-) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitleY("Events/[20 MeV/c^2]");
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
    	hi_ekpkpkmprot_MM2 = new H1F("hi_ekpkpkmprot_MM2", "hi_ekpkpkmprot_MM2", 70, -0.2, 0.5);
		hi_ekpkpkmprot_MM2.setTitle("MM2");
		hi_ekpkpkmprot_MM2.setTitleX("MM2(eK^+K^+K^-p) [GeV^2/c^4]");
		hi_ekpkpkmprot_MM2.setTitleY("Events");
		hi_ekpkpkmprot_MM2.setFillColor(LIGHTGREEN);
		
		hi_efkpckp_MM_req_km = new H1F("hi_efkpckp_MM_req_km", "hi_efkpckp_MM_req_km", 50, 1.6, 2.1);
		hi_efkpckp_MM_req_km.setTitle("MM requiring an additional K^- anywhere");
		hi_efkpckp_MM_req_km.setTitleX("MM(efk^+ck^+) (GeV)");
		hi_efkpckp_MM_req_km.setTitleY("count");
		hi_efkpckp_MM_req_km.setFillColor(LIGHTGREEN);
		
		hi_efkpckp_MM_req_ckm = new H1F("hi_efkpckp_MM_req_ckm", "hi_efkpckp_MM_req_ckm", 50, 1.6, 2.1);
		hi_efkpckp_MM_req_ckm.setTitle("MM requiring an additional K^- in CD");
		hi_efkpckp_MM_req_ckm.setTitleX("MM(efk^+ck^+) (GeV)");
		hi_efkpckp_MM_req_ckm.setTitleY("count");
		hi_efkpckp_MM_req_ckm.setFillColor(LIGHTGREEN);

		hi_ekpkp_MM_req_pim = new H1F("hi_ekpkp_MM_req_pim", "hi_ekpkp_MM_req_pim", 50, 1.6, 2.1);
		hi_ekpkp_MM_req_pim.setTitle("MM requiring an additional #pi^-");
		hi_ekpkp_MM_req_pim.setTitleX("MM(ek^+k^+) (GeV)");
		hi_ekpkp_MM_req_pim.setTitleY("count");
		hi_ekpkp_MM_req_pim.setFillColor(LIGHTGREEN);
		
		hi_efkpckp_MM_req_pim = new H1F("hi_efkpckp_MM_req_pim", "hi_efkpckp_MM_req_pim", 50, 1.6, 2.1);
		hi_efkpckp_MM_req_pim.setTitle("MM requiring an additional #pi^- anywhere");
		hi_efkpckp_MM_req_pim.setTitleX("MM(efk^+ck^+) (GeV)");
		hi_efkpckp_MM_req_pim.setTitleY("count");
		hi_efkpckp_MM_req_pim.setFillColor(LIGHTGREEN);
		
		hi_efkpckp_MM_req_cpim = new H1F("hi_efkpckp_MM_req_cpim", "hi_efkpckp_MM_req_cpim", 50, 1.6, 2.1);
		hi_efkpckp_MM_req_cpim.setTitle("MM requiring an additional #pi^- in CD");
		hi_efkpckp_MM_req_cpim.setTitleX("MM(efk^+ck^+) (GeV)");
		hi_efkpckp_MM_req_cpim.setTitleY("count");
		hi_efkpckp_MM_req_cpim.setFillColor(LIGHTGREEN);
		
		//FD particle vertex
		hi_pim_vz = new H1F("hi_pim_vz", "hi_pim_vz", 100, -20, 20);
		hi_pim_vz.setFillColor(LIGHTGREEN);
		hi_pip_vz = new H1F("hi_pip_vz", "hi_pip_vz", 100, -20, 20);
		hi_pip_vz.setFillColor(LIGHTGREEN);
		hi_kp_vz = new H1F("hi_kp_vz", "hi_kp_vz", 100, -20, 20);
		hi_kp_vz.setFillColor(LIGHTGREEN);
		hi_km_vz = new H1F("hi_km_vz", "hi_km_vz", 100, -20, 20);
		hi_km_vz.setFillColor(LIGHTGREEN);
		hi_prot_vz = new H1F("hi_prot_vz", "hi_prot_vz", 100, -20, 20);
		hi_prot_vz.setFillColor(LIGHTGREEN);
		
		//FD particle path
		hi_prot_FTOF1b_path = new H1F("hi_prot_FTOF1b_path", "hi_prot_FTOF1b_path", 200, 400, 900);
		hi_prot_FTOF1b_path.setFillColor(LIGHTGREEN);
		hi_pip_FTOF1b_path = new H1F("hi_pip_FTOF1b_path", "hi_pip_FTOF1b_path", 200, 400, 900);
		hi_pip_FTOF1b_path.setFillColor(LIGHTGREEN);
		hi_pim_FTOF1b_path = new H1F("hi_pim_FTOF1b_path", "hi_pim_FTOF1b_path", 200, 400, 900);
		hi_pim_FTOF1b_path.setFillColor(LIGHTGREEN);
		hi_kp_FTOF1b_path = new H1F("hi_kp_FTOF1b_path", "hi_kp_FTOF1b_path", 200, 400, 900);
		hi_kp_FTOF1b_path.setFillColor(LIGHTGREEN);
		hi_km_FTOF1b_path = new H1F("hi_km_FTOF1b_path", "hi_km_FTOF1b_path", 200, 400, 900);
		hi_km_FTOF1b_path.setFillColor(LIGHTGREEN);
		//FD particle time
		hi_prot_FTOF1b_t = new H1F("hi_prot_FTOF1b_t", "hi_prot_FTOF1b_t", 200, 0, 80);
		hi_prot_FTOF1b_t.setFillColor(LIGHTGREEN);
		hi_pip_FTOF1b_t = new H1F("hi_pip_FTOF1b_t", "hi_pip_FTOF1b_t", 200, 0, 80);
		hi_pip_FTOF1b_t.setFillColor(LIGHTGREEN);
		hi_pim_FTOF1b_t = new H1F("hi_pim_FTOF1b_t", "hi_pim_FTOF1b_t", 200, 0, 80);
		hi_pim_FTOF1b_t.setFillColor(LIGHTGREEN);
		hi_kp_FTOF1b_t = new H1F("hi_kp_FTOF1b_t", "hi_kp_FTOF1b_t", 200, 0, 80);
		hi_kp_FTOF1b_t.setFillColor(LIGHTGREEN);
		hi_km_FTOF1b_t = new H1F("hi_km_FTOF1b_t", "hi_km_FTOF1b_t", 200, 0, 80);
		hi_km_FTOF1b_t.setFillColor(LIGHTGREEN);
		
		// FTOF particles vetrex time vs momentum
		hi_pip_vt_p = new H2F("hi-pip-vt-p","hi-pip-vt-p", 100, 0, 8, 50, -0.5, 0.5);
		hi_pip_vt_p.setTitle("pip vt vs mom");
		hi_pip_vt_p.setTitleY("vt (ns)");
		hi_pip_vt_p.setTitleX("p (GeV)");
		hi_pim_vt_p = new H2F("hi_pim_vt_p","hi_pim_vt_p", 100, 0, 8, 50, -0.5, 0.5);
		hi_pim_vt_p.setTitle("pim vt vs mom");
		hi_pim_vt_p.setTitleY("vt (ns)");
		hi_pim_vt_p.setTitleX("p (GeV)");
		hi_kp_vt_p = new H2F("hi_kp_vt_p","hi_kp_vt_p", 100, 0, 8, 50, -0.5, 0.5);
		hi_kp_vt_p.setTitle("kp vt vs mom");
		hi_kp_vt_p.setTitleY("vt (ns)");
		hi_kp_vt_p.setTitleX("p (GeV)");
		hi_km_vt_p = new H2F("hi_km_vt_p","hi_km_vt_p", 100, 0, 8, 50, -0.5, 0.5);
		hi_km_vt_p.setTitle("km vt vs mom");
		hi_km_vt_p.setTitleY("vt (ns)");
		hi_km_vt_p.setTitleX("p (GeV)");
		hi_prot_vt_p = new H2F("hi_prot_vt_p","hi_prot_vt_p", 100, 0, 8, 50, -0.5, 0.5);
		hi_prot_vt_p.setTitle("prot vt vs mom");
		hi_prot_vt_p.setTitleY("vt (ns)");
		hi_prot_vt_p.setTitleX("p (GeV)");
		
		//CD particle vertex
		hi_pimc_vz = new H1F("hi_pimc_vz", "hi_pimc_vz", 100, -20, 20);
		hi_pimc_vz.setFillColor(LIGHTGREEN);
		hi_pipc_vz = new H1F("hi_pipc_vz", "hi_pipc_vz", 100, -20, 20);
		hi_pipc_vz.setFillColor(LIGHTGREEN);
		hi_kpc_vz = new H1F("hi_kpc_vz", "hi_kpc_vz", 100, -20, 20);
		hi_kpc_vz.setFillColor(LIGHTGREEN);
		hi_kmc_vz = new H1F("hi_km_vz", "hi_kmc_vz", 100, -20, 20);
		hi_kmc_vz.setFillColor(LIGHTGREEN);
		hi_protc_vz = new H1F("hi_protc_vz", "hi_protc_vz", 100, -20, 20);
		hi_protc_vz.setFillColor(LIGHTGREEN);
		
		//CD particle path
		hi_prot_CTOF_path = new H1F("hi_prot_CTOF_path", "hi_prot_CTOF_path", 200, 0, 100);
		hi_prot_CTOF_path.setFillColor(LIGHTGREEN);
		hi_pip_CTOF_path = new H1F("hi_pip_CTOF_path", "hi_pip_CTOF_path", 200, 0, 100);
		hi_pip_CTOF_path.setFillColor(LIGHTGREEN);
		hi_pim_CTOF_path = new H1F("hi_pim_CTOF_path", "hi_pim_CTOF_path", 200, 0, 100);
		hi_pim_CTOF_path.setFillColor(LIGHTGREEN);
		hi_kp_CTOF_path = new H1F("hi_kp_CTOF_path", "hi_kp_CTOF_path", 200, 0, 100);
		hi_kp_CTOF_path.setFillColor(LIGHTGREEN);
		hi_km_CTOF_path = new H1F("hi_km_CTOF_path", "hi_km_CTOF_path", 200, 0, 100);
		hi_km_CTOF_path.setFillColor(LIGHTGREEN);
		
		//CD particle time
		hi_prot_CTOF_t = new H1F("hi_prot_CTOF_t", "hi_prot_CTOF_t", 200, -5, 20);
		hi_prot_CTOF_t.setFillColor(LIGHTGREEN);
		hi_pip_CTOF_t = new H1F("hi_pip_CTOF_t", "hi_pip_CTOF_t", 200, -5, 20);
		hi_pip_CTOF_t.setFillColor(LIGHTGREEN);
		hi_pim_CTOF_t = new H1F("hi_pim_CTOF_t", "hi_pim_CTOF_t", 200, -5, 20);
		hi_pim_CTOF_t.setFillColor(LIGHTGREEN);
		hi_kp_CTOF_t = new H1F("hi_kp_CTOF_t", "hi_kp_CTOF_t", 200, -5, 20);
		hi_kp_CTOF_t.setFillColor(LIGHTGREEN);
		hi_km_CTOF_t = new H1F("hi_km_CTOF_t", "hi_km_CTOF_t", 200, -5, 20);
		hi_km_CTOF_t.setFillColor(LIGHTGREEN);		
		
		
		//CTOF particles vetrex time vs momentum
		hi_pipc_vt_p = new H2F("hi_pipc_vt_p","hi_pipc_vt_p", 50, 0, 3, 50, -0.5, 0.5);
		hi_pipc_vt_p.setTitle("pipc vt vs mom");
		hi_pipc_vt_p.setTitleY("vt (ns)");
		hi_pipc_vt_p.setTitleX("p (GeV)");
		
		hi_pimc_vt_p = new H2F("hi_pimc_vt_p","hi_pimc_vt_p", 50, 0, 3, 50, -0.5, 0.5);
		hi_pimc_vt_p.setTitle("pimc vt vs mom");
		hi_pimc_vt_p.setTitleY("vt (ns)");
		hi_pimc_vt_p.setTitleX("p (GeV)");
		
		hi_kpc_vt_p = new H2F("hi_kpc_vt_p","hi_kpc_vt_p", 50, 0, 3, 50, -0.5, 0.5);
		hi_kpc_vt_p.setTitle("kpc vt vs mom");
		hi_kpc_vt_p.setTitleY("vt (ns)");
		hi_kpc_vt_p.setTitleX("p (GeV)");
		
		hi_kmc_vt_p = new H2F("hi_kmc_vt_p","hi_kmc_vt_p", 50, 0, 3, 50, -0.5, 0.5);
		hi_kmc_vt_p.setTitle("kmc vt vs mom");
		hi_kmc_vt_p.setTitleY("vt (ns)");
		hi_kmc_vt_p.setTitleX("p (GeV)");
		
		hi_protc_vt_p = new H2F("hi_protc_vt_p","hi_protc_vt_p", 50, 0, 3, 50, -0.5, 0.5);
		hi_protc_vt_p.setTitle("protc vt vs mom");
		hi_protc_vt_p.setTitleY("vt (ns)");
		hi_protc_vt_p.setTitleX("p (GeV)");
		
		
		// FD particle beta vs momentum by charge
		//hi_FD_pos_beta_mom = new H2F("hi_FD_pos_beta_mom", "hi_FD_pos_beta_mom", 80, 0, 8, 50, 0.3, 1.1);
		hi_FD_pos_beta_mom = new H2F("H-FD-pos-beta-mom", "H-FD-pos-beta-mom", 250, 0.5, 5, 200, 0.85, 1.05); //80, 0.4, 5, 50, 0.85, 1.05
		hi_FD_pos_beta_mom.setTitle("POS  #beta vs mom");
		hi_FD_pos_beta_mom.setTitleX("p (GeV)");
		hi_FD_pos_beta_mom.setTitleY("FTB #beta");
		//hi_FD_neg_beta_mom = new H2F("hi_FD_neg_beta_mom", "hi_FD_neg_beta_mom", 80, 0, 8, 50, 0.3, 1.1);
		hi_FD_neg_beta_mom = new H2F("H-FD-neg-beta-mom", "H-FD-neg-beta-mom", 250, 0.5, 5, 200, 0.85, 1.05); //80, 0.4, 5, 50, 0.85, 1.05
		hi_FD_neg_beta_mom.setTitle("NEG  #beta vs mom");
		hi_FD_neg_beta_mom.setTitleX("p (GeV)");
		hi_FD_neg_beta_mom.setTitleY("FTB #beta");
		hi_FD_neutral_beta_mom = new H2F("hi_FD_neutral_beta_mom", "hi_FD_neutral_beta_mom", 80, 0, 8, 50, 0.3, 1.1);
		hi_FD_neutral_beta_mom.setTitle("NEUTRAL  #beta vs mom");
		hi_FD_neutral_beta_mom.setTitleX("p (GeV)");
		hi_FD_neutral_beta_mom.setTitleY("FTB #beta");
		hi_fd_pos_mass = new H1F("hi_fd_pos_mass", "hi_fd_pos_mass", 800, -0.5, 4.5);
		hi_fd_pos_mass.setFillColor(LIGHTGREEN);
		hi_FD_pos_mass_mom = new H2F("hi_FD_pos_mass_mom", "hi_FD_pos_mass_mom", 100, 0, 7, 800, -0.5, 4.5);
		hi_FD_pos_mass_mom.setTitle("POS Mass^2 vs mom");
		hi_FD_pos_mass_mom.setTitleX("p (GeV)");
		hi_FD_pos_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_fd_neg_mass = new H1F("hi_fd_neg_mass", "hi_fd_neg_mass", 500, -0.5, 2);
		hi_fd_neg_mass.setFillColor(LIGHTGREEN);
		hi_FD_neg_mass_mom = new H2F("hi_FD_neg_mass_mom", "hi_FD_neg_mass_mom", 100, 0, 7, 500, -0.5, 2);
		hi_FD_neg_mass_mom.setTitle("NEG Mass^2 vs mom");
		hi_FD_neg_mass_mom.setTitleX("p (GeV)");
		hi_FD_neg_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_FD_neutral_mass_mom = new H2F("hi_FD_neutral_mass_mom", "hi_FD_neutral_mass_mom", 100, 0, 7, 500, -0.5, 2);
		hi_FD_neutral_mass_mom.setTitle("NEUTRAL Mass^2 vs mom");
		hi_FD_neutral_mass_mom.setTitleX("p (GeV)");
		hi_FD_neutral_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_FD_pos_mass_the = new H2F("hi_FD_pos_mass_the", "hi_FD_pos_mass_the", 100, 0, 80, 800, -0.5, 4.5);
		hi_FD_pos_mass_the.setTitle("POS Mass^2 vs #theta");
		hi_FD_pos_mass_the.setTitleX("#theta (^o)");
		hi_FD_pos_mass_the.setTitleY("M^2 (GeV^2)");
		hi_FD_neg_mass_the = new H2F("hi_FD_neg_mass_the", "hi_FD_neg_mass_the", 100, 0, 80, 500, -0.5, 2);
		hi_FD_neg_mass_the.setTitle("NEG Mass^2 vs #theta");
		hi_FD_neg_mass_the.setTitleX("#theta (^o)");
		hi_FD_neg_mass_the.setTitleY("M^2 (GeV^2)");
		hi_FD_neutral_mass_the = new H2F("hi_FD_neutral_mass_the", "hi_FD_neutral_mass_the", 100, 0, 80, 500, -0.5, 2);
		hi_FD_neutral_mass_the.setTitle("NEUTRAL Mass^2 vs #theta");
		hi_FD_neutral_mass_the.setTitleX("#theta (^o)");
		hi_FD_neutral_mass_the.setTitleY("M^2 (GeV^2)");
		hi_FD_pos_mass_phi = new H2F("hi_FD_pos_mass_phi", "hi_FD_pos_mass_phi", 100, -180, 180, 800, -0.5, 4.5);
		hi_FD_pos_mass_phi.setTitle("POS Mass^2 vs #phi");
		hi_FD_pos_mass_phi.setTitleX("#phi (^o)");
		hi_FD_pos_mass_phi.setTitleY("M^2 (GeV^2)");
		hi_FD_neg_mass_phi = new H2F("hi_FD_neg_mass_phi", "hi_FD_neg_mass_phi", 100, -180, 180, 500, -0.5, 2);
		hi_FD_neg_mass_phi.setTitle("NEG Mass^2 vs #phi");
		hi_FD_neg_mass_phi.setTitleX("#phi (^o)");
		hi_FD_neg_mass_phi.setTitleY("M^2 (GeV^2)");
		hi_FD_neutral_mass_phi = new H2F("hi_FD_neutral_mass_phi", "hi_FD_neutral_mass_phi", 100, -180, 180, 500, -0.5, 2);
		hi_FD_neutral_mass_phi.setTitle("NEUTRAL Mass^2 vs #phi");
		hi_FD_neutral_mass_phi.setTitleX("#phi (^o)");
		hi_FD_neutral_mass_phi.setTitleY("M^2 (GeV^2)");
		
		
		// CD particle beta vs mom by charge
		hi_CD_pos_beta_mom = new H2F("hi_CD_pos_beta_mom", "hi_CD_pos_beta_mom", 50, 0, 2.0, 50, 0.1, 1.1);
		hi_CD_pos_beta_mom.setTitle("POS  #beta vs mom");
		hi_CD_pos_beta_mom.setTitleX("p (GeV)");
		hi_CD_pos_beta_mom.setTitleY("FTB #beta");
		hi_CD_neg_beta_mom = new H2F("hi_CD_neg_beta_mom", "hi_CD_neg_beta_mom", 50, 0, 2.0, 50, 0.1, 1.1);
		hi_CD_neg_beta_mom.setTitle("NEG  #beta vs mom");
		hi_CD_neg_beta_mom.setTitleX("p (GeV)");
		hi_CD_neg_beta_mom.setTitleY("FTB #beta");
		hi_CD_neutral_beta_mom = new H2F("hi_CD_neutral_beta_mom", "hi_CD_neutral_beta_mom", 50, 0, 2.0, 50, 0.1, 1.1);
		hi_CD_neutral_beta_mom.setTitle("NEUTRAL  #beta vs mom");
		hi_CD_neutral_beta_mom.setTitleX("p (GeV)");
		hi_CD_neutral_beta_mom.setTitleY("FTB #beta");
		hi_cd_pos_mass = new H1F("hi_cd_pos_mass", "hi_cd_pos_mass", 800, -0.5, 4.5);
		hi_cd_pos_mass.setFillColor(LIGHTGREEN);
		hi_CD_pos_mass_mom = new H2F("hi_CD_pos_mass_mom", "hi_CD_pos_mass_mom", 50, 0, 2.0, 800, -0.5, 4.5);
		hi_CD_pos_mass_mom.setTitle("POS Mass^2 vs mom");
		hi_CD_pos_mass_mom.setTitleX("p (GeV)");
		hi_CD_pos_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_cd_neg_mass = new H1F("hi_fd_neg_mass", "hi_fd_neg_mass", 500, -0.5, 2.0);
		hi_cd_neg_mass.setFillColor(LIGHTGREEN);
		hi_CD_neg_mass_mom = new H2F("hi_CD_neg_mass_mom", "hi_CD_neg_mass_mom", 50, 0, 2.0, 500, -0.5, 2.0);
		hi_CD_neg_mass_mom.setTitle("NEG Mass^2 vs mom");
		hi_CD_neg_mass_mom.setTitleX("p (GeV)");
		hi_CD_neg_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_CD_neutral_mass_mom = new H2F("hi_CD_neutral_mass_mom", "hi_CD_neutral_mass_mom", 50, 0, 2.0, 500, -0.5, 2.0);
		hi_CD_neutral_mass_mom.setTitle("NEUTRAL Mass^2 vs mom");
		hi_CD_neutral_mass_mom.setTitleX("p (GeV)");
		hi_CD_neutral_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_CD_pos_mass_the = new H2F("hi_CD_pos_mass_the", "hi_CD_pos_mass_the", 100, 30, 130, 800, -0.5, 4.5);
		hi_CD_pos_mass_the.setTitle("POS Mass^2 vs #theta");
		hi_CD_pos_mass_the.setTitleX("#theta (^o)");
		hi_CD_pos_mass_the.setTitleY("M^2 (GeV^2)");
		hi_CD_neg_mass_the = new H2F("hi_CD_neg_mass_the", "hi_CD_neg_mass_the", 100, 30, 130, 500, -0.5, 2);
		hi_CD_neg_mass_the.setTitle("NEG Mass^2 vs #theta");
		hi_CD_neg_mass_the.setTitleX("#theta (^o)");
		hi_CD_neg_mass_the.setTitleY("M^2 (GeV^2)");
		hi_CD_neutral_mass_the = new H2F("hi_CD_neutral_mass_the", "hi_CD_neutral_mass_the", 100, 30, 130, 500, -0.5, 2);
		hi_CD_neutral_mass_the.setTitle("NEUTRAL Mass^2 vs #theta");
		hi_CD_neutral_mass_the.setTitleX("#theta (^o)");
		hi_CD_neutral_mass_the.setTitleY("M^2 (GeV^2)");
		hi_CD_pos_mass_phi = new H2F("hi_CD_pos_mass_phi", "hi_CD_pos_mass_phi", 150, -180, 180, 800, -0.5, 4.5);
		hi_CD_pos_mass_phi.setTitle("POS Mass^2 vs #phi");
		hi_CD_pos_mass_phi.setTitleX("#phi (^o)");
		hi_CD_pos_mass_phi.setTitleY("M^2 (GeV^2)");
		hi_CD_neg_mass_phi = new H2F("hi_CD_neg_mass_phi", "hi_CD_neg_mass_phi", 150, -180, 180, 500, -0.5, 2);
		hi_CD_neg_mass_phi.setTitle("NEG Mass^2 vs #phi");
		hi_CD_neg_mass_phi.setTitleX("#phi (^o)");
		hi_CD_neg_mass_phi.setTitleY("M^2 (GeV^2)");
		hi_CD_neutral_mass_phi = new H2F("hi_CD_neutral_mass_phi", "hi_CD_neutral_mass_phi", 150, -180, 180, 500, -0.5, 2);
		hi_CD_neutral_mass_phi.setTitle("NEUTRAL Mass^2 vs #phi");
		hi_CD_neutral_mass_phi.setTitleX("#phi (^o)");
		hi_CD_neutral_mass_phi.setTitleY("M^2 (GeV^2)");
		

		// rec electron
		dg_rec_electron = new DataGroup(3, 5);
		dg_rec_electron.addDataSet(hi_FT_e_p_the, 0);
		dg_rec_electron.addDataSet(hi_FT_e_p_f, 1);
		dg_rec_electron.addDataSet(hi_FT_e_t_f, 2);
		dg_rec_electron.addDataSet(hi_rec_e_dp, 3);
		dg_rec_electron.addDataSet(fn_rec_e_dp, 3);
		dg_rec_electron.addDataSet(hi_rec_e_dtheta, 4);
		dg_rec_electron.addDataSet(fn_rec_e_dtheta, 4);
		dg_rec_electron.addDataSet(hi_rec_e_dphi, 5);
		dg_rec_electron.addDataSet(fn_rec_e_dphi, 5);
		dg_rec_electron.addDataSet(hi_FT_Q2, 6);
		dg_rec_electron.addDataSet(hi_FT_W, 7);
		dg_rec_electron.addDataSet(hi_virphoton, 8);
		dg_rec_electron.addDataSet(hi_FT_W_Q2, 9);
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
		dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkp_nocorr, 2);
		dg_rec_xi.addDataSet(fn_xi_no_corr, 2);
		dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkpkm_nocorr, 3);
		//dg_rec_xi.addDataSet(f1_sigma, 1);
		dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkp, 4);
		dg_rec_xi.addDataSet(f1_xi, 4);
		dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkpkm, 5);
		//dg_rec_xi.addDataSet(hi_ekpkpkm_MM_ekpkpkm_nocorr, 3);
		//dg_rec_xi.addDataSet(f1_lambda, 3);
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
		dg_mm.addDataSet(hi_ekpkp_MM_req_pim, 0);
		dg_mm.addDataSet(hi_efkpckp_MM_req_pim, 1);
		dg_mm.addDataSet(hi_efkpckp_MM_req_cpim, 2);
		dg_mm.addDataSet(hi_ekpkpkm_MM_ekpkp, 3);
		dg_mm.addDataSet(hi_efkpckp_MM_req_km, 4);
		dg_mm.addDataSet(hi_efkpckp_MM_req_ckm, 5);
		dg_mm.addDataSet(hi_ekpkpkmprot_MM2, 6);

		dg_mm_ekpkpkm = new DataGroup(3, 2);
		dg_mm_ekpkpkm.addDataSet(hi_efkpfkpfkm_MM_efkpfkp, 0);
		dg_mm_ekpkpkm.addDataSet(hi_efkpckpfkm_MM_efkpckp, 1);
		dg_mm_ekpkpkm.addDataSet(hi_eckpckpfkm_MM_eckpckp, 2);
		dg_mm_ekpkpkm.addDataSet(hi_efkpfkpckm_MM_efkpfkp, 3);
		dg_mm_ekpkpkm.addDataSet(hi_efkpckpckm_MM_efkpckp, 4);
		dg_mm_ekpkpkm.addDataSet(hi_eckpckpckm_MM_eckpckp, 5);

		dg_mm_ekpkpkm1 = new DataGroup(3, 2);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpfkpfkm_MM_efkpfkpfkm, 0);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpckpfkm_MM_efkpckpfkm_nocorn, 1);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpckpfkm_MM_efkpckpfkm, 1);
		dg_mm_ekpkpkm1.addDataSet(hi_eckpckpfkm_MM_eckpckpfkm, 2);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpfkpckm_MM_efkpfkpckm, 3);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpckpckm_MM_efkpckpckm, 4);
		dg_mm_ekpkpkm1.addDataSet(hi_eckpckpckm_MM_eckpckpckm, 5);

		dg_mm_scatter = new DataGroup(3, 2);
		dg_mm_scatter.addDataSet(hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm, 0);
		dg_mm_scatter.addDataSet(hi_efkpckpfkm_MM_efkpckp_MM_efkpckpfkm, 1);
		dg_mm_scatter.addDataSet(hi_eckpckpfkm_MM_eckpckp_MM_eckpckpfkm, 2);
		dg_mm_scatter.addDataSet(hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm, 3);
		dg_mm_scatter.addDataSet(hi_efkpckpckm_MM_efkpckp_MM_efkpckpckm, 4);
		dg_mm_scatter.addDataSet(hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm, 5);

		dg_mm_ekpkpkm_lam_evnt = new DataGroup(3, 2);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_efkpfkpfkm_MM_efkpfkp_lam_evnt, 0);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_efkpckpfkm_MM_efkpckp_lam_evnt, 1);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_eckpckpfkm_MM_eckpckp_lam_evnt, 2);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_efkpfkpckm_MM_efkpfkp_lam_evnt, 3);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_efkpckpckm_MM_efkpckp_lam_evnt, 4);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_eckpckpckm_MM_eckpckp_lam_evnt, 5);

		dg_m_kmlambda = new DataGroup(3, 2);
		dg_m_kmlambda.addDataSet(hi_efkpfkpfkm_IM_kmlambda, 0);
		dg_m_kmlambda.addDataSet(hi_efkpckpfkm_IM_kmlambda, 1);
		dg_m_kmlambda.addDataSet(hi_eckpckpfkm_IM_kmlambda, 2);
		dg_m_kmlambda.addDataSet(hi_efkpfkpckm_IM_kmlambda, 3);
		dg_m_kmlambda.addDataSet(hi_efkpckpckm_IM_kmlambda, 4);
		dg_m_kmlambda.addDataSet(hi_eckpckpckm_IM_kmlambda, 5);


		dg_fkpckp_pthphi = new DataGroup(4,4);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltap, 0);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltatheta, 1);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltaphi, 2);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltap_p, 4);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltatheta_p, 5);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltaphi_p, 6);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltap_withdpdtcut, 0);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltap_theta, 8);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltatheta_theta, 9);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltaphi_theta, 10);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltatheta_withdpdtcut, 1);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltap_phi, 12);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltatheta_phi, 13);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltaphi_phi, 14);
		dg_fkpckp_pthphi.addDataSet(hi_fkpckp_deltaphi_withdpdtcut, 2);
		
		
		
		// verticx time vs p

		F1D fn_fkp_vtime_ucut = new F1D("fn-fkp-vtime-ucut", "0.05+exp(-1.0*x+0.4/x)", 1.2, 8);
		F1D fn_fkp_vtime_lcut = new F1D("fn-fkp-vtime-lcut", "-(0.05+exp(-1.0*x+0.4/x))", 1.2, 8);
		F1D fn_ckp_vtime_ucut = new F1D("fn-ckp-vtime-ucut", "0.05+exp(-6.0*x+0.15/x)", 0.4, 3);
		F1D fn_ckp_vtime_lcut = new F1D("fn-ckp-vtime-lcut", "-(0.05+exp(-3.5*x+0.15/x))", 0.4, 3);
		F1D fn_fkm_vtime_ucut = new F1D("fn-fkm-vtime-ucut", "0.05+exp(-1.0*x+0.4/x)", 1.2, 8);
		F1D fn_fkm_vtime_lcut = new F1D("fn-fkm-vtime-lcut", "-(0.05+exp(-1.0*x+0.4/x))", 1.2, 8);
		F1D fn_ckm_vtime_ucut = new F1D("fn-ckm-vtime-ucut", "0.05+exp(-6.0*x+0.15/x)", 0.4, 3);
		F1D fn_ckm_vtime_lcut = new F1D("fn-ckm-vtime-lcut", "-(0.05+exp(-3.5*x+0.15/x))", 0.4, 3);

		dg_vtime = new DataGroup(5, 2);
		dg_vtime.addDataSet(hi_pip_vt_p, 0);
		dg_vtime.addDataSet(hi_pim_vt_p, 1);
		dg_vtime.addDataSet(hi_kp_vt_p, 2);
		dg_vtime.addDataSet(fn_fkp_vtime_ucut, 2);
		dg_vtime.addDataSet(fn_fkp_vtime_lcut, 2);
		dg_vtime.addDataSet(hi_km_vt_p, 3);
		dg_vtime.addDataSet(fn_fkm_vtime_ucut, 3);
		dg_vtime.addDataSet(fn_fkm_vtime_lcut, 3);
		dg_vtime.addDataSet(hi_prot_vt_p, 4);
		dg_vtime.addDataSet(hi_pipc_vt_p, 5);
		dg_vtime.addDataSet(hi_pimc_vt_p, 6);
		dg_vtime.addDataSet(hi_kpc_vt_p, 7);
		dg_vtime.addDataSet(fn_ckp_vtime_ucut, 7);
		dg_vtime.addDataSet(fn_ckp_vtime_lcut, 7);
		dg_vtime.addDataSet(hi_kmc_vt_p, 8);
		dg_vtime.addDataSet(fn_ckm_vtime_ucut, 8);
		dg_vtime.addDataSet(fn_ckm_vtime_lcut, 8);
		dg_vtime.addDataSet(hi_protc_vt_p, 9);
		
		dg_fdtof = new DataGroup(5, 2);
		dg_fdtof.addDataSet(hi_FD_pos_beta_mom, 0);
		dg_fdtof.addDataSet(F_prot_beta_mom, 0);
		dg_fdtof.addDataSet(F_pip_beta_mom,0);
		dg_fdtof.addDataSet(F_kp_beta_mom,0);
		dg_fdtof.addDataSet(hi_FD_pos_mass_mom, 1);
		dg_fdtof.addDataSet(hi_FD_pos_mass_the, 2);
		dg_fdtof.addDataSet(hi_FD_pos_mass_phi, 3);
		dg_fdtof.addDataSet(hi_fd_pos_mass, 4);
		dg_fdtof.addDataSet(hi_FD_neg_beta_mom, 5);
		dg_fdtof.addDataSet(F_prot_beta_mom, 5);
		dg_fdtof.addDataSet(F_pip_beta_mom, 5);
		dg_fdtof.addDataSet(F_kp_beta_mom, 5);
		dg_fdtof.addDataSet(hi_FD_neg_mass_mom, 6);
		dg_fdtof.addDataSet(hi_FD_neg_mass_the, 7);
		dg_fdtof.addDataSet(hi_FD_neg_mass_phi, 8);
		dg_fdtof.addDataSet(hi_fd_neg_mass, 9);

		dg_cdtof = new DataGroup(5, 2);
		dg_cdtof.addDataSet(hi_CD_pos_beta_mom, 0);
		dg_cdtof.addDataSet(F_prot_beta_mom, 0);
		dg_cdtof.addDataSet(F_pip_beta_mom,0);
		dg_cdtof.addDataSet(F_kp_beta_mom,0);
		dg_cdtof.addDataSet(hi_CD_pos_mass_mom, 1);
		dg_cdtof.addDataSet(hi_CD_pos_mass_the, 2);
		dg_cdtof.addDataSet(hi_CD_pos_mass_phi, 3);
		dg_cdtof.addDataSet(hi_cd_pos_mass, 4);
		dg_cdtof.addDataSet(hi_CD_neg_beta_mom, 5);
		dg_cdtof.addDataSet(F_prot_beta_mom, 5);
		dg_cdtof.addDataSet(F_pip_beta_mom, 5);
		dg_cdtof.addDataSet(F_kp_beta_mom, 5);
		dg_cdtof.addDataSet(hi_CD_neg_mass_mom, 6);
		dg_cdtof.addDataSet(hi_CD_neg_mass_the, 7);
		dg_cdtof.addDataSet(hi_CD_neg_mass_phi, 8);
		dg_cdtof.addDataSet(hi_cd_neg_mass, 9);

		// CD/FD hadron particles vz
		dg_vz = new DataGroup(5, 2);
		dg_vz.addDataSet(hi_pip_vz, 0);
		dg_vz.addDataSet(hi_pim_vz, 1);
		dg_vz.addDataSet(hi_kp_vz, 2);
		dg_vz.addDataSet(hi_km_vz, 3);
		dg_vz.addDataSet(hi_prot_vz, 4);
		dg_vz.addDataSet(hi_pipc_vz, 5);
		dg_vz.addDataSet(hi_pimc_vz, 6);
		dg_vz.addDataSet(hi_kpc_vz, 7);
		dg_vz.addDataSet(hi_kmc_vz, 8);
		dg_vz.addDataSet(hi_protc_vz, 9);

		//TOF_t, 

		dg_tof_t = new DataGroup(5, 2);
		dg_tof_t.addDataSet(hi_pip_FTOF1b_t, 0);
		dg_tof_t.addDataSet(hi_pim_FTOF1b_t, 1);
		dg_tof_t.addDataSet(hi_kp_FTOF1b_t, 2);
		dg_tof_t.addDataSet(hi_km_FTOF1b_t, 3);
		dg_tof_t.addDataSet(hi_prot_FTOF1b_t, 4);
		dg_tof_t.addDataSet(hi_pip_CTOF_t, 5);
		dg_tof_t.addDataSet(hi_pim_CTOF_t, 6);
		dg_tof_t.addDataSet(hi_kp_CTOF_t, 7);
		dg_tof_t.addDataSet(hi_km_CTOF_t, 8);
		dg_tof_t.addDataSet(hi_prot_CTOF_t, 9);

		//TOF_path
		dg_tof_path = new DataGroup(5, 2);
		dg_tof_path.addDataSet(hi_pip_FTOF1b_path, 0);
		dg_tof_path.addDataSet(hi_pim_FTOF1b_path, 1);
		dg_tof_path.addDataSet(hi_kp_FTOF1b_path, 2);
		dg_tof_path.addDataSet(hi_km_FTOF1b_path, 3);
		dg_tof_path.addDataSet(hi_prot_FTOF1b_path, 4);
		dg_tof_path.addDataSet(hi_pip_CTOF_path, 5);
		dg_tof_path.addDataSet(hi_pim_CTOF_path, 6);
		dg_tof_path.addDataSet(hi_kp_CTOF_path, 7);
		dg_tof_path.addDataSet(hi_km_CTOF_path, 8);
		dg_tof_path.addDataSet(hi_prot_CTOF_path, 9);

		//CD protons and pips
		hi_cd_prot_p = new H1F("hi_cd_prot_p", "hi-cd-prot-p", 50, 0, 2.0);
		hi_cd_prot_p.setTitleX("p (GeV)");
		hi_cd_prot_p.setTitleY("count");
		hi_cd_prot_p.setFillColor(LIGHTGREEN);
		hi_cd_pip_p = new H1F("hi_cd_prot_p", "hi-cd-pip-p", 50, 0, 2.0);
		hi_cd_pip_p.setTitleX("p (GeV)");
		hi_cd_pip_p.setTitleY("count");
		hi_cd_pip_p.setFillColor(LIGHTGREEN);
		hi_cd_prot_theta = new H1F("hi_cd_prot_theta", "hi-cd-prot-theta", 100, 30, 130);
		hi_cd_prot_theta.setTitleX("#theta (^o)");
		hi_cd_prot_theta.setTitleY("count");
		hi_cd_prot_theta.setFillColor(LIGHTGREEN);
		hi_cd_pip_theta = new H1F("hi_cd_prot_theta", "hi-cd-pip-theta", 100, 30, 130);
		hi_cd_pip_theta.setTitleX("#theta (^o)");
		hi_cd_pip_theta.setTitleY("count");
		hi_cd_pip_theta.setFillColor(LIGHTGREEN);
		hi_cd_prot_p_theta = new H2F("hi-cd-prot-p-theta", "hi-cd-prot-p-theta", 100, 30, 130, 50, 0, 2.0);
		hi_cd_prot_p_theta.setTitleX("#theta (^o)");
		hi_cd_prot_p_theta.setTitleY("p (GeV)");
		hi_cd_pip_p_theta = new H2F("hi-cd-pip-p-theta", "hi-cd-pip-p-theta", 100, 30, 130, 50, 0, 2.0);
		hi_cd_pip_p_theta.setTitleX("#theta (^o)");
		hi_cd_pip_p_theta.setTitleY("p (GeV)");

		hi_W_cd_pro_the = new H2F("hi-W-cd-pro-the", "hi-W-cd-pro-the", 100, 3, 5, 100, 30, 130);//100, 3, 5
		hi_W_cd_pro_the.setTitleX("W (GeV)");
		hi_W_cd_pro_the.setTitleY("CD protons #theta (^o)");

		dg_cdPart = new DataGroup(3,3);
		dg_cdPart.addDataSet(hi_cd_prot_theta, 0);
		dg_cdPart.addDataSet(hi_cd_prot_p_theta, 1);
		dg_cdPart.addDataSet(hi_cd_prot_p, 2);
		dg_cdPart.addDataSet(hi_cd_pip_theta, 3);
		dg_cdPart.addDataSet(hi_cd_pip_p_theta, 4);
		dg_cdPart.addDataSet(hi_cd_pip_p, 5);

		dg_cdPart.addDataSet(hi_W_cd_pro_the, 6);

		
	} // end of ft_ana()

	public void fillMCPartBank(DataBank mcbank, DataBank genconfig){

		int mcrunnum = genconfig.getInt('run',0);
        int mcevnum  =  genconfig.getInt('event',0);

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

		// labeling reconstructed kps with momentum
		Particle Vslowkp;
		Particle Vfastkp;

		if(mckps.get(0).p() > mckps.get(1).p()){
			Vfastkp = new Particle(mckps.get(0));
			Vslowkp = new Particle(mckps.get(1));
		} else {
			Vfastkp = new Particle(mckps.get(1));
			Vslowkp = new Particle(mckps.get(0));
		}


		LorentzVector lv_q = new LorentzVector(VB); lv_q.sub(mcels.get(0).vector());


		//for mass calculation
		LorentzVector lv_fkpskp = new LorentzVector();lv_fkpskp.add(Vfastkp.vector());lv_fkpskp.add(Vslowkp.vector());
		double m_fkpskp = lv_fkpskp.mass();
		LorentzVector lv_fkpskpkm = new LorentzVector();lv_fkpskpkm.add(lv_fkpskp);lv_fkpskpkm.add(mckms.get(0).vector());
		double m_fkpskpkm = lv_fkpskpkm.mass();
		LorentzVector lv_fkpkm = new LorentzVector();lv_fkpkm.add(Vfastkp.vector());lv_fkpkm.add(mckms.get(0).vector());
		double m_fkpkm = lv_fkpkm.mass();
		LorentzVector lv_skpkm = new LorentzVector();lv_skpkm.add(Vslowkp.vector());lv_fkpkm.add(mckms.get(0).vector());
		double m_skpkm = lv_skpkm.mass();
		//missing mass calculation
		LorentzVector lv_mm_e = new LorentzVector();lv_mm_e.add(lv_q);lv_mm_e.add(VT);
		double mm_e = lv_mm_e.mass();
		LorentzVector lv_mm_fkp = new LorentzVector();lv_mm_fkp.add(lv_mm_e);lv_mm_fkp.sub(Vfastkp.vector());
		double mm_fkp = lv_mm_fkp.mass();
		LorentzVector lv_mm_skp = new LorentzVector();lv_mm_skp.add(lv_mm_e);lv_mm_skp.sub(Vslowkp.vector());
		double mm_skp = lv_mm_skp.mass();
		LorentzVector lv_mm_km = new LorentzVector();lv_mm_km.add(lv_mm_e);lv_mm_km.sub(mckms.get(0).vector());
		double mm_km =  lv_mm_km.mass();
		LorentzVector lv_mm_fkpskp = new LorentzVector();lv_mm_fkpskp.add(lv_mm_fkp);lv_mm_fkpskp.sub(Vslowkp.vector());
		double mm_fkpskp = lv_mm_fkpskp.mass();
		LorentzVector lv_mm_fkpkm = new LorentzVector();lv_mm_fkpkm.add(lv_mm_fkp);lv_mm_fkpkm.sub(mckms.get(0).vector());
		double mm_fkpkm = lv_mm_fkpkm.mass();
		LorentzVector lv_mm_skpkm = new LorentzVector();lv_mm_skpkm.add(lv_mm_skp);lv_mm_skpkm.sub(mckms.get(0).vector());
		double mm_skpkm = lv_mm_skpkm.mass();
		LorentzVector lv_mm_fkpskpkm = new LorentzVector();lv_mm_fkpskpkm.add(lv_mm_fkpskp);lv_mm_fkpskpkm.sub(mckms.get(0).vector());
		double mm_fkpskpkm = lv_mm_fkpskpkm.mass();

		//electron kinematics
		double q2 = -lv_q.mass2();
		double nu = VB.e()-mcels.get(0).e();
		double x  = q2 / (2 * PDGDatabase.getParticleById(2212).mass() * nu);
		double W  = Math.pow(Math.pow(PDGDatabase.getParticleById(2212).mass(),2)+2*PDGDatabase.getParticleById(2212).mass()*nu - q2, 0.5);

		// lab e kinematics
		double beam_e = VB.e();
		double e_px = mcels.get(0).px();
		double e_py = mcels.get(0).py();
		double e_pz = mcels.get(0).pz();
		double e_p = mcels.get(0).p();
		double e_e = mcels.get(0).e();
		double e_vx = mcels.get(0).vx();
		double e_vy = mcels.get(0).vy();
		double e_vz = mcels.get(0).vz();
		double e_theta = Math.toDegrees(mcels.get(0).theta());
		double e_phi = Math.toDegrees(mcels.get(0).phi());
		// lab fkp kinematics
		double fkp_px = Vfastkp.px();
		double fkp_py = Vfastkp.py();
	    double fkp_pz = Vfastkp.pz();
		double fkp_p = Vfastkp.p();
		double fkp_e = Vfastkp.e();
		double fkp_vx = Vfastkp.vx();
		double fkp_vy = Vfastkp.vy();
		double fkp_vz = Vfastkp.vz();
		double fkp_theta = Math.toDegrees(Vfastkp.theta());
		double fkp_phi = Math.toDegrees(Vfastkp.phi());
		// lab skp kinematics
		double skp_px = Vslowkp.px();
		double skp_py = Vslowkp.py();
		double skp_pz = Vslowkp.pz();
		double skp_p = Vslowkp.p();
		double skp_e = Vslowkp.e();
		double skp_vx = Vslowkp.vx();
		double skp_vy = Vslowkp.vy();
		double skp_vz = Vslowkp.vz();
		double skp_theta = Math.toDegrees(Vslowkp.theta());
		double skp_phi = Math.toDegrees(Vslowkp.phi());
		// lab km kinematics
		double km_px = mckms.get(0).vector().px();
		double km_py = mckms.get(0).vector().py();
		double km_pz = mckms.get(0).vector().pz();
		double km_p = mckms.get(0).vector().p();
		double km_e = mckms.get(0).vector().e();
		double km_vx = mckms.get(0).vx();
		double km_vy = mckms.get(0).vy();
		double km_vz = mckms.get(0).vz();
		double km_theta = Math.toDegrees(mckms.get(0).vector().theta());
		double km_phi = Math.toDegrees(mckms.get(0).vector().phi());

		// append event to next line of the text file
		//file.append(runnum+" "+evnum+" "+helicity+" ");
		mcoutFile.append(mcrunnum+" "+mcevnum+" "+beam_e+" "+q2+" "+nu+" "+x+" "+W+" ");
		mcoutFile.append(e_px+" "+e_py+" "+e_pz+" "+e_p+" "+e_e+" "+e_vx+" "+e_vy+" "+e_vz+" "+e_theta+" "+e_phi+" ");
		mcoutFile.append(fkp_px+" "+fkp_py+" "+fkp_pz+" "+fkp_p+" "+fkp_e+" "+fkp_vx+" "+fkp_vy+" "+fkp_vz+" "+fkp_theta+" "+fkp_phi+" ");
		mcoutFile.append(skp_px+" "+skp_py+" "+skp_pz+" "+skp_p+" "+skp_e+" "+skp_vx+" "+skp_vy+" "+skp_vz+" "+skp_theta+" "+skp_phi+" ");
		mcoutFile.append(km_px+" "+km_py+" "+km_pz+" "+km_p+" "+km_e+" "+km_vx+" "+km_vy+" "+km_vz+" "+km_theta+" "+km_phi+" ");
		mcoutFile.append(m_fkpskp+" "+m_fkpkm+" "+m_skpkm+" "+m_fkpskpkm+" ");
		mcoutFile.append(mm_e+" "+mm_fkp+" "+mm_skp+" "+mm_km+" "+mm_fkpskp+" "+mm_fkpkm+" "+mm_skpkm+" "+mm_fkpskpkm+"\n");
		//outFile.append(mm_fkpaspipskpkm+" "+mm_fkpskpaspipkm+" "+mm_fkpskpkmaspim+" "+mm_fkpaspipskpaspipkm+" "+mm_fkpaspipskpkmaspim+" "+mm_fkpskpaspipkmaspim+" "+mm_fkpaspipskpaspipkmaspim+"\n");
		//outFile.append("\n");
		//println(); println();
		//print("1:runnum, 2:evnum, 3:beam_e, 4:q2, 5:nu, 6:x, 7:W, 8:e_px, 9:e_py, 10:e_pz, 11:e_p, 12:e_e, 13:e_vx, 14:e_vy, 15:e_vz, 16:e_theta, 17:e_phi,");
		//print("14:fkp_px, 15:fkp_py, 16:fkp_pz, 17:fkp_p, 18:fkp_e, 19:fkp_vx, 20:fkp_vy, 21:fkp_vz, 22:fkp_theta, 23:fkp_phi,");
		//print("24:skp_px, 25:skp_py, 26:skp_pz, 27:skp_p, 28:skp_e, 29:skp_vx, 30:skp_vy, 31:skp_vz, 32:skp_theta, 33:skp_phi,");
		//print("34:km_px, 35:km_py, 36:km_pz, 37:km_p, 38:km_e, 39:km_vx, 40:km_vy, 41:km_vz, 42:km_theta, 43:km_phi,");


		
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

//FD smearing functions (From Giovanni/Florian) to match MC and data resolution derived using elastic event
	public LorentzVector smear(LorentzVector Vector){

	  	Random r = new Random();
	   	double inM = Vector.mass();
	   	double sP  = Vector.p();
	   	double sTh = Vector.theta();
	    double sPh = Vector.phi();
	   	double sThD = Math.toDegrees(sTh);
	   	double momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
	    double momS2 = 0.02*sThD;
	    double momR  = 0.01 * Math.sqrt(  Math.pow(momS1*sP, 2) + Math.pow(momS2, 2));
	    momR *= 2.0; // <- to match data resolution

	    double theS1 = 0.004*sThD + 0.1;
	    double theS2 = 0;
	    double theR  = Math.sqrt(Math.pow(theS1*Math.sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + Math.pow(theS2,2));
	    theR *= 2.5; // <- to match data resolution

	    double phiS1 = 0.85-0.015*sThD;
	    double phiS2 = 0.17-0.003*sThD;
	    double phiR  = Math.sqrt(Math.pow(phiS1*Math.sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + Math.pow(phiS2,2) );
	    phiR *= 3.5; // <- to match data resolution

	    sPh +=Math.toRadians(phiR)  * r.nextGaussian();
	    sTh += Math.toRadians(theR)  * r.nextGaussian();
	    sP  += momR  * r.nextGaussian() * Vector.p() ; 
	    LorentzVector FinalVector = new LorentzVector();
	    double px = sP* Math.sin(sTh)* Math.cos(sPh);
		double py = sP* Math.sin(sTh)* Math.sin(sPh);
		double pz = sP*Math.cos(sTh);
		FinalVector.setPxPyPzE(px, py, pz, Math.sqrt( sP*sP + inM*inM ) );
		return FinalVector;
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
			short status = (short) Math.abs(recFT.getShort("status", k));
			boolean inFT = (partstatus >= 1000 && partstatus < 2000 && recftstatus < 0.0 && (int) (status / 1000) == 1);
			e_mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			e_the = (float) Math.toDegrees(Math.acos(pz / e_mom));
			if (ftbpid == 11 && q == -1 && inFT  && e_mom > 0.1 && e_mom < 4.5 && e_the < 4.5 && e_the > 2.5) { //&& e_mom > 0.5 && e_mom < 4.5
				NFTElec++;
				Velectron = new Particle(ftbpid, px, py, pz, vx, vy, vz);
				// applying FT electron energy correction
				Vftel = electron_energy_correction(Velectron);

				ftels.add(Vftel);

				if (runType == "mc"){
					//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
					Random r = new Random();
    				double esmearFactor = 2.5*r.nextGaussian()/100;
    				double smearedP = Velectron.p()*(1+esmearFactor);
					Velectron.setP(smearedP);
				}

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
				// correction from /volatile/clas12/users/devita/2pi/rga-pass1v0/skim3/
				//double eftNew = eft-0.03689+0.1412*eft-0.04316*eft*eft+0.007046*eft*eft*eft-0.0004055*eft*eft*eft*eft;
				//latest correction from /volatile/clas12/users/devita/ft/vertex/dst/recon/005038/
				double eftNew = eft+0.1574-0.01178*eft-0.007345*eft*eft+0.002909*eft*eft*eft-0.0002243*eft*eft*eft*eft;
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
				
				//Ve = new LorentzVector(px, py, pz, eft);
				Ve = Vftel.vector();
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
	}

	
	public void makeOthers(DataBank recbank,  DataBank recFTbank, DataBank recSCBank) {
		int nkp = 0;
		int npip = 0;
		int nfdpip = 0;
		int nfdwithcutpip = 0;
		int ncdpip = 0;
		int ncdwithcutpip = 0;
		int nfdandcdpip = 0;
		
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

		// momentum smearing factor for MC reconstruction
		//Random rand = new Random();
    	//double smearFactor = rand.nextGaussian()/100;
		
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
			float vt = recFTbank.getFloat("vt", k);
			if (status<0) status = -status;
			boolean inDC = (status >= 2000 && status < 4000);
			boolean inCD = (status >= 4000);
			float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			float the = (float) Math.toDegrees(Math.acos(pz / mom));
			float phi = (float) Math.toDegrees(Math.atan2(py, px));
			//float FTBmass = mom * mom * (1 / ( ftbbe * ftbbe ) - 1);
			float FTBmass = mom * mom * ((1 -  ftbbe * ftbbe)/ ( ftbbe * ftbbe ));
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
							//pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_beta * 29.98f) - STT - pip_vz/ (pip_beta * 29.98f);
							//pip_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							//float pip_FTOF1b_TOF = (float) pip_FTOF1b_t - STT - pip_vz/ (pip_beta * 29.98f);

							//(STT + pip_vz/ (ftbbe* 29.98f))

							//float pip_FTOF1b_TOF = (float) pip_FTOF1b_t - STT - pip_vz/ (ftbbe* 29.98f); 
							float pip_FTOF1b_TOF = (float) pip_FTOF1b_t - recFTbank.getFloat("vt", k);

							//difference in the measured and computed vertex times 
							pip_FTOF1b_vt = (float) pip_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass())));
							if (Math.abs(pip_FTOF1b_vt) < 0.6) { //0.4//0.5
								nfdwithcutpip++;
								//Vpip = new LorentzVector(pip_px*pipCor, pip_py*pipCor, pip_pz*pipCor, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								//Vpip = new LorentzVector(pip_px, pip_py, pip_pz, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								Vpip = new Particle(ftbpid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);
								fpips.add(Vpip);
								pips.add(Vpip);
								//hi_pip_vt_p.fill(pip_FTOF1b_vt, pip_mom);

								//assign the PID=211 and -211 events the kaon mass.
								Particle Vpipaskp = new Particle();
								Vpipaskp = new Particle(321 , pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);

								fkps.add(Vpipaskp);
								kps.add(Vpipaskp);
							}
							hi_pip_vt_p.fill(pip_mom, pip_FTOF1b_vt);
							hi_pip_vz.fill(pip_vz);
							hi_pip_FTOF1b_t.fill(pip_FTOF1b_TOF);
							hi_pip_FTOF1b_path.fill(pip_FTOF1b_path);
							
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
							//pip_CTOF_vt = pip_CTOF_t - pip_CTOF_path / (pipc_beta * 29.98f) - STT - pip_vz/ (pipc_beta * 29.98f);
							//pip_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							//float pip_CTOF_TOF = (float) pip_CTOF_t - STT - pip_vz/ (pipc_beta * 29.98f);

							//float pip_CTOF_TOF = (float) pip_CTOF_t - STT - pip_vz/ (ftbbe * 29.98f);
							float pip_CTOF_TOF = (float) pip_CTOF_t - recFTbank.getFloat("vt", k);
							//difference in the measured and computed vertex times 
							pip_CTOF_vt = (float) pip_CTOF_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass())));
							if (Math.abs(pip_CTOF_vt) < 0.5) {//0.4
								ncdwithcutpip++;
								//Vpipc = new LorentzVector(pip_px*pipCor, pip_py*pipCor, pip_pz*pipCor, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								
								Particle cdpart = new Particle(ftbpid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);
								Vpipc = new Particle(ftbpid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);	
								//Vpipc = cd_momentum_correction(cdpart);
								pips.add(Vpipc);
								cpips.add(Vpipc);
								//hi_pipc_vt_p.fill(pip_FTOF1b_vt, pip_mom);

								//assign the PID=211 and -211 events the kaon mass.
								Particle Vpipaskp = new Particle();
								Vpipaskp = new Particle(321 , pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);

								ckps.add(Vpipaskp);
								kps.add(Vpipaskp);
							}
							hi_pipc_vt_p.fill(pip_mom, pip_CTOF_vt);
							hi_pipc_vz.fill(pip_vz);
							hi_pip_CTOF_t.fill(pip_CTOF_t - STT - pip_vz/ (pipc_beta * 29.98f));
							hi_pip_CTOF_path.fill(pip_CTOF_path);
							
							
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
							//pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_beta * 29.98f) - STT - pim_vz/ (pim_beta * 29.98f);
							//pim_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							//float pim_FTOF1b_TOF = (float) pim_FTOF1b_t - STT - pim_vz/ (pim_beta * 29.98f);

							//float pim_FTOF1b_TOF = (float) pim_FTOF1b_t - STT - pim_vz/ (ftbbe * 29.98f);
							float pim_FTOF1b_TOF = (float) pim_FTOF1b_t - recFTbank.getFloat("vt", k);

							//difference in the measured and computed vertex times 
							pim_FTOF1b_vt = (float) pim_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass())));
							if (Math.abs(pim_FTOF1b_vt) < 0.6) {//0.5
								//Vpim = new LorentzVector(pim_px*pimCor, pim_py*pimCor, pim_pz*pimCor, Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass()));
								Vpim = new Particle(ftbpid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);

								if (runType == "mc"){
									//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
									//Random rand = new Random();
    								//double smearFactor = rand.nextGaussian()/100;
    								double smearedP = Vpim.p()*(1+smearFactor);
									Vpim.setP(smearedP);
								}

								pims.add(Vpim);
								fpims.add(Vpim);
								//hi_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);

								//assign the PID=211 and -211 events the kaon mass.
								Particle Vpimaskm = new Particle();
								Vpimaskm = new Particle(-321 , pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);

								fkms.add(Vpimaskm);
								kms.add(Vpimaskm);
							}
							hi_pim_vt_p.fill(pim_mom, pim_FTOF1b_vt);
							hi_pim_vz.fill(pim_vz);
							hi_pim_FTOF1b_t.fill(pim_FTOF1b_t - STT - pim_vz/ (pim_beta * 29.98f));
							hi_pim_FTOF1b_path.fill(pim_FTOF1b_path);
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
							//float pim_CTOF_TOF = (float) pim_CTOF_t - STT - pim_vz/ (pimc_beta * 29.98f);

							//float pim_CTOF_TOF = (float) pim_CTOF_t - STT - pim_vz/ (ftbbe * 29.98f);
							float pim_CTOF_TOF = (float) pim_CTOF_t - recFTbank.getFloat("vt", k);

							//difference in the measured and computed vertex times 
							pim_CTOF_vt = (float) pim_CTOF_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass())));
							if (Math.abs(pim_CTOF_vt) < 0.5) {//0.4
								//Vpimc = new LorentzVector(pim_px*pimCor, pim_py*pimCor, pim_pz*pimCor, Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass()));
								Particle cdpart = new Particle(ftbpid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);
								//Vpimc = cd_momentum_correction(cdpart);
								Vpimc = new Particle(ftbpid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);		
								pims.add(Vpimc);
								cpims.add(Vpimc);
								//hi_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);

								//assign the PID=211 and -211 events the kaon mass.
								Particle Vpimaskm = new Particle();
								Vpimaskm = new Particle(-321 , pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);

								ckms.add(Vpimaskm);
								kms.add(Vpimaskm);
							}
							hi_pimc_vt_p.fill(pim_mom, pim_CTOF_vt);
							hi_pimc_vz.fill(pim_vz);
							hi_pim_CTOF_t.fill(pim_CTOF_t - STT - pim_vz/ (pimc_beta * 29.98f));
							hi_pim_CTOF_path.fill(pim_CTOF_path);
							
							
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
							//kp_FTOF1b_vt = kp_FTOF1b_t - kp_FTOF1b_path / (kp_beta * 29.98f) - STT - kp_vz/ (kp_beta * 29.98f);
							//kp_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							//float kp_FTOF1b_TOF = (float) kp_FTOF1b_t - STT - kp_vz/ (kp_beta * 29.98f);

							//float kp_FTOF1b_TOF = (float) kp_FTOF1b_t - STT - kp_vz/ (ftbbe * 29.98f);
							float kp_FTOF1b_TOF = (float) kp_FTOF1b_t - recFTbank.getFloat("vt", k);

							//difference in the measured and computed vertex times 
							kp_FTOF1b_vt = (float) kp_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass())));
							boolean fdkpvzCut = vz > -10 && vz < 2;
							float fkp_vtime_cut = Math.abs(0.05+Math.exp(-1.0*mom+0.4/mom));
							boolean fdkpCuts = (Math.abs(kp_FTOF1b_vt) <  fkp_vtime_cut && kp_FTOF1b_TOF > 20 && kp_FTOF1b_TOF < 55 && ftbbe > 0.4 && ftbbe < 1.05 && mom > 0.4 && mom < Eb && fdkpvzCut);
							if (fdkpCuts ) {//fdkpCuts // && kp_mom < 2.8 //fdkpvzCut && 
								Vkp = new Particle(ftbpid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);

								if (runType == "mc"){
									//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
									//Random rand = new Random();
    								//double smearFactor = rand.nextGaussian()/100;
    								double smearedP = Vkp.p()*(1+smearFactor);
									Vkp.setP(smearedP);
								}

								//fkps.add(Vkp);
								//kps.add(Vkp);
								//hi_kp_vt_p.fill(kp_FTOF1b_vt, kp_mom);
							}
							hi_kp_vt_p.fill(kp_mom, kp_FTOF1b_vt);
							hi_kp_vz.fill(kp_vz);
							hi_kp_FTOF1b_t.fill(kp_FTOF1b_t - STT - kp_vz/ (kp_beta * 29.98f));
							hi_kp_FTOF1b_path.fill(kp_FTOF1b_path);
						} // kp from FTOF panal 1b
					}
					
					if (recSCBank.getByte("detector", r) == 4 && inCD  && cdChi2pidCut ) { //&& cdkpMass2Cut
						
						if (recSCBank.getShort("pindex", r) == kp_part_ind ) {
							kp_CTOF_pad = recSCBank.getShort("component", r);
							kp_CTOF_t = recSCBank.getFloat("time", r);
							kp_CTOF_path = recSCBank.getFloat("path", r);
							float kpc_beta = kp_mom / (float) Math.sqrt(kp_mom * kp_mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass());
							//kp_CTOF_vt = kp_CTOF_t - kp_CTOF_path / (kpc_beta * 29.98f) - STT - kp_vz/ (kpc_beta * 29.98f);
							//kp_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							boolean cdkpvzCut = vz > -10 && vz < 2;
							//float kp_CTOF_TOF = (float) kp_CTOF_t - STT - kp_vz/ (kpc_beta * 29.98f);

							//float kp_CTOF_TOF = (float) kp_CTOF_t - STT - kp_vz/ (ftbbe * 29.98f);
							float kp_CTOF_TOF = (float) kp_CTOF_t - recFTbank.getFloat("vt", k);

							//difference in the measured and computed vertex times 
							kp_CTOF_vt = (float) kp_CTOF_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass())));
							float ckp_vtime_ucut = 0.05+Math.exp(-6.0*mom+0.15/mom);
							float ckp_vtime_lcut = -(0.05+Math.exp(-3.5*mom+0.15/mom));
							boolean cdkpCut = kp_CTOF_vt < ckp_vtime_ucut && kp_CTOF_vt > ckp_vtime_lcut && kp_CTOF_TOF > 0.5 && kp_CTOF_TOF < 4 && ftbbe > 0.2 && ftbbe < 1.05 && mom > 0.2 && mom < 3.0 && cdkpvzCut;
							if (cdkpCut ) {// cdkpCut //cdkpvzCut && Math.abs(kp_CTOF_vt) < 0.4
								 /*
								Particle cdpart = new Particle(ftbpid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);
								Vkpc = cd_momentum_correction(cdpart); // */
								Vkpc = new Particle(ftbpid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);
								//kps.add(Vkpc);
								//ckps.add(Vkpc);
								//hi_kp_vt_p.fill(kp_FTOF1b_vt, kp_mom);
							}
							hi_kpc_vt_p.fill(kp_mom, kp_CTOF_vt);
							hi_kpc_vz.fill(kp_vz);
							hi_kp_CTOF_t.fill(kp_CTOF_t - STT - kp_vz/ (kpc_beta * 29.98f));
							hi_kp_CTOF_path.fill(kp_CTOF_path);
							
							
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
							//km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_beta * 29.98f) - STT - km_vz/ (km_beta * 29.98f);
							//km_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							boolean fdkmvzCut = vz > -10 && vz < 2;

							//float km_FTOF1b_TOF = (float) km_FTOF1b_t - STT - km_vz/ (km_beta * 29.98f);

							//float km_FTOF1b_TOF = (float) km_FTOF1b_t - STT - km_vz/ (ftbbe * 29.98f);
							float km_FTOF1b_TOF = (float) km_FTOF1b_t - recFTbank.getFloat("vt", k);

							boolean km_FTOF1b_TOFCut = 22 < km_FTOF1b_TOF && km_FTOF1b_TOF < 28;
							//difference in the measured and computed vertex times 
							km_FTOF1b_vt = (float) km_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass())));
							float fkm_vtime_cut = Math.abs(0.05+Math.exp(-1.0*mom+0.4/mom));
							boolean fdkmCut = Math.abs(km_FTOF1b_vt) < fkm_vtime_cut && km_FTOF1b_TOF > 20 && km_FTOF1b_TOF < 35 && ftbbe > 0.4 && ftbbe < 1.05 && mom > 0.4 && mom < Eb && fdkmvzCut;
							if (fdkmCut) {//fdkmCut //&& fdkmvzCut && km_FTOF1b_TOFCut Math.abs(km_FTOF1b_vt) < 0.5
								Vkm = new Particle(ftbpid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);

								if (runType == "mc"){
									//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
									//Random rand = new Random();
    								//double smearFactor = rand.nextGaussian()/100;
    								double smearedP = Vkm.p()*(1+smearFactor);
									Vkm.setP(smearedP);
								}

								//fkms.add(Vkm);
								//kms.add(Vkm);
								//hi_km_vt_p.fill(km_FTOF1b_vt, km_mom);
							}
							hi_km_vt_p.fill(km_mom, km_FTOF1b_vt);
							hi_km_vz.fill(km_vz);
							hi_km_FTOF1b_t.fill(km_FTOF1b_TOF);
							hi_km_FTOF1b_path.fill(km_FTOF1b_path);
						} // km from FTOF panal 1b
					}
					
					if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut ) {
						
						if (recSCBank.getShort("pindex", r) == km_part_ind ) {
							km_CTOF_pad = recSCBank.getShort("component", r);
							km_CTOF_t = recSCBank.getFloat("time", r);
							km_CTOF_path = recSCBank.getFloat("path", r);
							float kmc_beta = km_mom / (float) Math.sqrt(km_mom * km_mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass());
							//km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_ftb_beta * 29.98f) - STT - km_vz/ (km_ftb_beta * 29.98f);
							//km_CTOF_vt = km_CTOF_t - km_CTOF_path / (kmc_beta * 29.98f) - STT - km_vz/ (kmc_beta * 29.98f);
							//km_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							boolean cdkmvzCut = vz > -10 && vz < 2;
							//float km_CTOF_TOF = (float) km_CTOF_t - STT - km_vz/ (kmc_beta * 29.98f);

							//float km_CTOF_TOF = (float) km_CTOF_t - STT - km_vz/ (ftbbe * 29.98f);
							float km_CTOF_TOF = (float) km_CTOF_t - recFTbank.getFloat("vt", k);

							//difference in the measured and computed vertex times 
							km_CTOF_vt = (float) km_CTOF_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass())));
							float ckm_vtime_ucut = 0.05+Math.exp(-6.0*mom+0.15/mom);
							float ckm_vtime_lcut = -(0.05+Math.exp(-3.5*mom+0.15/mom));
							boolean cdkmCut = km_CTOF_vt < ckm_vtime_ucut && km_CTOF_vt > ckm_vtime_lcut && km_CTOF_TOF > 0.5 && km_CTOF_TOF < 4 && ftbbe > 0.2 && ftbbe < 1.05 && mom > 0.2 && mom < 3.0 && cdkmvzCut;
							if (cdkmCut) { //&& cdkmvzCut Math.abs(km_CTOF_vt) < 0.4 
								
								Particle cdpart = new Particle(ftbpid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);
								//Vkmc = cd_momentum_correction(cdpart);
								Vkmc = new Particle(ftbpid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);		
								//kms.add(Vkmc);
								//ckms.add(Vkmc);
								//hi_km_vt_p.fill(km_FTOF1b_vt, km_mom);
							}
							hi_kmc_vt_p.fill(km_mom, km_CTOF_vt);
							hi_kmc_vz.fill(km_vz);
							hi_km_CTOF_t.fill(km_CTOF_t - STT - km_vz/ (kmc_beta * 29.98f));
							hi_km_CTOF_path.fill(km_CTOF_path);
							
							
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
							//prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_beta * 29.98f) - STT - prot_vz/ (prot_beta * 29.98f);
							//prot_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							//float prot_FTOF1b_TOF = (float) prot_FTOF1b_t - STT - prot_vz/ (prot_beta * 29.98f);

							//float prot_FTOF1b_TOF = (float) prot_FTOF1b_t - STT - prot_vz/ (ftbbe * 29.98f);
							float prot_FTOF1b_TOF = (float) prot_FTOF1b_t - recFTbank.getFloat("vt", k);

							boolean prot_FTOF1b_TOFCut = 22 < prot_FTOF1b_TOF && prot_FTOF1b_TOF < 32;
							//difference in the measured and computed vertex times 
							prot_FTOF1b_vt = (float) prot_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass())));
							if (Math.abs(prot_FTOF1b_vt) < 0.6 ) {// && prot_FTOF1b_TOFCut
								//Vprot = new LorentzVector(prot_px * proCor, prot_py* proCor, prot_pz* proCor, Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass()));
								Vprot = new Particle(ftbpid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);
								fprots.add(Vprot);
								prots.add(Vprot);
								//hi_prot_vt_p.fill(prot_FTOF1b_vt, prot_mom);
							}
							hi_prot_vt_p.fill(prot_mom, prot_FTOF1b_vt);
							hi_prot_vz.fill(prot_vz);
							hi_prot_FTOF1b_t.fill(prot_FTOF1b_t - STT - prot_vz/ (prot_beta * 29.98f));
							hi_prot_FTOF1b_path.fill(prot_FTOF1b_path);
						} // prot from FTOF panal 1b
					}
					
					if (recSCBank.getByte("detector", r) == 4 && inCD && cdChi2pidCut ) {
						
						if (recSCBank.getShort("pindex", r) == prot_part_ind ) {
							prot_CTOF_pad = recSCBank.getShort("component", r);
							prot_CTOF_t = recSCBank.getFloat("time", r);
							prot_CTOF_path = recSCBank.getFloat("path", r);
							float protc_beta = prot_mom / (float) Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass());
							//prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_ftb_beta * 29.98f) - STT - prot_vz/ (prot_ftb_beta * 29.98f);
							//prot_CTOF_vt = prot_CTOF_t - prot_CTOF_path / (protc_beta * 29.98f) - STT - prot_vz/ (protc_beta * 29.98f);
							//prot_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							//float prot_CTOF_TOF = (float) prot_CTOF_t - STT - prot_vz/ (protc_beta * 29.98f);

							//float prot_CTOF_TOF = (float) prot_CTOF_t - STT - prot_vz/ (ftbbe * 29.98f);
							float prot_CTOF_TOF = (float) prot_CTOF_t - recFTbank.getFloat("vt", k);

							//difference in the measured and computed vertex times 
							prot_CTOF_vt = (float) prot_CTOF_TOF * (1 - Math.sqrt((mom*mom + FTBmass)/(mom*mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass())));
							if (Math.abs(prot_CTOF_vt) < 0.5) {
								//Vprotc = new LorentzVector(prot_px * proCor, prot_py* proCor, prot_pz* proCor, Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass()));
								 /*
								Particle cdpart = new Particle(ftbpid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);
								Vprotc = cd_momentum_correction(cdpart); //*/
								Vprotc = new Particle(ftbpid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);		
								cprots.add(Vprotc);
								prots.add(Vprotc);
								//hi_prot_vt_p.fill(prot_FTOF1b_vt, prot_mom);
							}
							hi_protc_vt_p.fill(prot_mom, prot_CTOF_vt);
							hi_protc_vz.fill(prot_vz);
							hi_prot_CTOF_t.fill(prot_CTOF_t - STT - prot_vz/ (protc_beta * 29.98f));
							hi_prot_CTOF_path.fill(prot_CTOF_path);
							
							
						}
						
						
					} //CTOF


					
					
				}
			
			}
			
			
			//fillFTOF(recbank, recSCBank);
			
			
			
			if ( q > 0 && ((ftbpid == 45 && vz > -9 && vz < 6) || (ftbpid == 2212 && vz > -9 && vz < 6) || (ftbpid == 321 && vz > -9 && vz < 6) || (ftbpid == 211 && vz > -9 && vz < 6) || (ftbpid == -11 && vz > -9 && vz < 6)) ) {
				if (inDC && fdChi2pidCut && mom > 0.6) {//&& fdChi2pidCut
					hi_fd_pos_mass.fill(FTBmass);
					hi_FD_pos_beta_mom.fill(mom, ftbbe);
					hi_FD_pos_mass_mom.fill(mom, FTBmass);
					hi_FD_pos_mass_the.fill(the, FTBmass);
					hi_FD_pos_mass_phi.fill(phi, FTBmass);		
				}
				if (inCD && cdChi2pidCut && mom > 0.2) {//&& cdChi2pidCut
					hi_cd_pos_mass.fill(FTBmass);
					hi_CD_pos_beta_mom.fill(mom, ftbbe);
					hi_CD_pos_mass_mom.fill(mom, FTBmass);
					hi_CD_pos_mass_the.fill(the, FTBmass);
					hi_CD_pos_mass_phi.fill(phi, FTBmass);	

					if(pid == 2212){
						hi_cd_prot_p.fill(mom);
						hi_cd_prot_theta.fill(the);
						hi_cd_prot_p_theta.fill(the, mom);

						hi_W_cd_pro_the.fill(e_W, the);
					}
					if(pid == 211){
						hi_cd_pip_p.fill(mom);
						hi_cd_pip_theta.fill(the);
						hi_cd_pip_p_theta.fill(the, mom);
					}				
				}
					
			}
			if ( q < 0 && ((ftbpid == -2212 && vz > -9 && vz < 6) || (ftbpid == -321 && vz > -9 && vz < 6) || (ftbpid == -211 && vz > -9 && vz < 6))) {
				if (inDC && fdChi2pidCut && mom > 0.6) {//&& fdChi2pidCut
					hi_fd_neg_mass.fill(FTBmass);
					hi_FD_neg_beta_mom.fill(mom, ftbbe);
					hi_FD_neg_mass_mom.fill(mom, FTBmass);
					hi_FD_neg_mass_the.fill(the, FTBmass);
					hi_FD_neg_mass_phi.fill(phi, FTBmass);		
				}
				if (inCD && cdChi2pidCut && mom > 0.2) {//&& cdChi2pidCut
					hi_cd_neg_mass.fill(FTBmass);
					hi_CD_neg_beta_mom.fill(mom, ftbbe);
					hi_CD_neg_mass_mom.fill(mom, FTBmass);
					hi_CD_neg_mass_the.fill(the, FTBmass);
					hi_CD_neg_mass_phi.fill(phi, FTBmass);				
				}

			}
			if (q == 0 ) {
				if (inDC ) {//&& fdChi2pidCut
					hi_FD_neutral_beta_mom.fill(mom, ftbbe);
					hi_FD_neutral_mass_mom.fill(mom, FTBmass);
					hi_FD_neutral_mass_the.fill(the, FTBmass);
					hi_FD_neg_mass_phi.fill(phi, FTBmass);		
				}
				if (inCD ) {//&& cdChi2pidCut
					hi_CD_neutral_beta_mom.fill(mom, ftbbe);
					hi_CD_neutral_mass_mom.fill(mom, FTBmass);
					hi_CD_neutral_mass_the.fill(the, FTBmass);
					hi_CD_neutral_mass_phi.fill(phi, FTBmass);				
				}
			}
			
		} // FOR LOOp
		
//		if(npip >= 2 && nfdpip == 2) {
//		System.out.println("found total pip in this event :: " + npip +" with " + nfdpip + " pips tracks in FD " + nfdwithcutpip +" with cut & " + ncdpip + " pip tracks in CD " + ncdwithcutpip + " with cut" );
//		}

	
	} //MAKEOTHER
	
	public boolean select_ekpkp_req_pim() {
		boolean res = false;
		/*
		if( found_eFT && kps.size() == 2 && pims.size() == 1) {
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(kps.get(0).vector());Vmissekpkp.sub(kps.get(1).vector());
			ekpkp_MM_req_pim = (float)Vmissekpkp.mass();
			res = true;
		}
		*/
		return res;
		
	}
	
	public boolean select_efkpckp_req_pim() {
		boolean res = false;
		/*
		if( found_eFT && fkps.size() >= 1  && ckps.size() >= 1 && pims.size() >= 1) {
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			efkpckp_MM_req_pim = (float)Vmissekpkp.mass();
			res = true;
		}
		*/
		return res;
	}
	
	public boolean select_efkpckp_req_cpim() {
		boolean res = false;
		/*
		if( found_eFT && fkps.size() >= 1  && ckps.size() >= 1 && cpims.size() >= 1) {
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			efkpckp_MM_req_cpim = (float)Vmissekpkp.mass();
			res = true;
		}
		*/
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
	public boolean select_efkpfkpfkm(){
		boolean res = false;
		efkpfkpfkm_found_lambda = false;
		if(found_eFT && fkps.size() == 2 && fkms.size() == 1 ){

			double avg_vz = (double)(fkps.get(0).vz() + fkps.get(1).vz() + fkms.get(0).vz())/3;
			//Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			Vftel_corrected = ftels.get(0);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Vftel_corrected.vector());Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(fkps.get(1).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			efkpfkpfkm_MM_efkpfkp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Vftel_corrected.vector());Vmissekpkpkm.sub(fkps.get(0).vector());Vmissekpkpkm.sub(fkps.get(1).vector());Vmissekpkpkm.sub(fkms.get(0).vector());
			efkpfkpfkm_MM_efkpfkpfkm = (float)Vmissekpkpkm.mass();

			if (Math.abs(efkpfkpfkm_MM_efkpfkpfkm - 1.132) <= 0.150 ){//Math.abs(efkpfkpfkm_MM_efkpfkpfkm - 1.126) <= 0.153 (3sigma cut 3*0.051)//efkpfkpfkm_MM_efkpfkpfkm < 1.17 && efkpfkpfkm_MM_efkpfkpfkm > 1.06(to selectLambdaK-decay only/excluds 2030)//efkpfkpfkm_MM_efkpfkpfkm < 1.28 && efkpfkpfkm_MM_efkpfkpfkm > 1.06(to include SigmaK- decay from 2030 as well)// efkpfkpfkm_MM_efkpfkpfkm < 1.27 && efkpfkpfkm_MM_efkpfkpfkm > 1.05 //Math.abs(efkpfkpfkm_MM_efkpfkpfkm - 1.115683) <= 0.052 
				efkpfkpfkm_found_lambda = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedlambda = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.115683f*1.115683f));
				LorentzVector VxifromLambda = new LorentzVector(0, 0, 0, 0);
				VxifromLambda.add(Vconstrainedlambda);VxifromLambda.add(fkms.get(0).vector());
				efkpfkpfkm_IM_kmlambda = (float)VxifromLambda.mass();
				
			}

			res = true;
		}
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


    public boolean select_efkpckpfkm(){
		boolean res = false;
		if(found_eFT && fkps.size() == 1&& ckps.size()==1 && fkms.size() == 1 ){


			double avg_vz0 = (double)(fkps.get(0).vz() + ckps.get(0).vz() + fkms.get(0).vz())/3;
			Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz0); //ftels.get(0).vector()
			LorentzVector lv_missekpkpkm_nocorn = new LorentzVector(0, 0, 0, 0);
			lv_missekpkpkm_nocorn.add(VT);lv_missekpkpkm_nocorn.add(VB);lv_missekpkpkm_nocorn.sub(Vftel_corrected.vector());lv_missekpkpkm_nocorn.sub(fkps.get(0).vector()); lv_missekpkpkm_nocorn.sub(ckps.get(0).vector());lv_missekpkpkm_nocorn.sub(fkms.get(0).vector());
			//float efkpckpfkm_MM_efkpckpfkm_nocorn;
			efkpckpfkm_MM_efkpckpfkm_nocorn = (float)lv_missekpkpkm_nocorn.mass();
			//hi_efkpckpfkm_MM_efkpckpfkm_nocorn.fill(efkpckpfkm_MM_efkpckpfkm_nocorn);
		   	/*
            boolean correctPhi = false;
            boolean correctTheta = false;
            boolean correctP = false;
            //*/
            //   /*
            boolean correctPhi = true;
            boolean correctTheta = true;
            boolean correctP = true;
            // */
            /*
            boolean correctPhi_ita2 = false;
            boolean correctTheta_ita2 = false;
            //*/

        //     /*
            boolean correctPhi_ita2 = true;
            boolean correctTheta_ita2 = true;
            //*/

            if(correctPhi ){ // correct CD Kp Phi

                ckps.get(0).vector().vect().setMagThetaPhi(ckps.get(0).p(), ckps.get(0).theta(), Math.toRadians(cdpart_phiCorrn(ckps.get(0))));
                
            }

            if(correctTheta){
                
                ckps.get(0).setTheta(Math.toRadians(cdpart_thetaCorrn(ckps.get(0))));
            }

            if(correctP){
                ckps.get(0).setP(cdpart_pCorrn(ckps.get(0)));
            }

            if(correctPhi_ita2){ 

                ckps.get(0).vector().vect().setMagThetaPhi(ckps.get(0).p(), ckps.get(0).theta(), Math.toRadians(cdpart_phiCorrn_it2(ckps.get(0))));
                
            }

            if(correctTheta_ita2){
                ckps.get(0).setTheta(Math.toRadians(cdpart_thetaCorrn_ita2(ckps.get(0))));
            }

			double avg_vz = (double)(fkps.get(0).vz() + ckps.get(0).vz() + fkms.get(0).vz())/3;
			Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);

			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Vftel_corrected.vector());Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			efkpckpfkm_MM_efkpckp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Vftel_corrected.vector());Vmissekpkpkm.sub(fkps.get(0).vector());Vmissekpkpkm.sub(ckps.get(0).vector());Vmissekpkpkm.sub(fkms.get(0).vector());
			efkpckpfkm_MM_efkpckpfkm = (float)Vmissekpkpkm.mass();

			if (Math.abs(efkpckpfkm_MM_efkpckpfkm - 1.115683) <= 0.052){//Math.abs(efkpfkpfkm_MM_efkpfkpfkm - 1.132) <= 0.150 //Math.abs(efkpckpfkm_MM_efkpckpfkm - 1.115683) <= 0.052
				efkpckpfkm_found_lambda = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedlambda = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.115683f*1.115683f));
				LorentzVector VxifromLambda = new LorentzVector(0, 0, 0, 0);
				VxifromLambda.add(Vconstrainedlambda);VxifromLambda.add(fkms.get(0).vector());
				efkpckpfkm_IM_kmlambda = (float)VxifromLambda.mass();
			}

			res = true;
		}
		return res;
	}

	public boolean select_eckpckpfkm(){
		boolean res = false;
		if(found_eFT && ckps.size()>=2 && fkms.size() >= 1 ){
			double avg_vz = (double)(ckps.get(0).vz() + ckps.get(1).vz() + fkms.get(0).vz())/3;
			Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Vftel_corrected.vector());Vmissekpkp.sub(ckps.get(0).vector());Vmissekpkp.sub(ckps.get(1).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			eckpckpfkm_MM_eckpckp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Vftel_corrected.vector());Vmissekpkpkm.sub(ckps.get(0).vector());Vmissekpkpkm.sub(ckps.get(1).vector());Vmissekpkpkm.sub(fkms.get(0).vector());
			eckpckpfkm_MM_eckpckpfkm = (float)Vmissekpkpkm.mass();

			if (Math.abs(eckpckpfkm_MM_eckpckpfkm - 1.115683) <= 0.052){
				eckpckpfkm_found_lambda = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedlambda = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.115683f*1.115683f));
				LorentzVector VxifromLambda = new LorentzVector(0, 0, 0, 0);
				VxifromLambda.add(Vconstrainedlambda);VxifromLambda.add(fkms.get(0).vector());
				eckpckpfkm_IM_kmlambda = (float)VxifromLambda.mass();
				
			}

			res = true;
		}
		return res;
	}

	public boolean select_efkpfkpckm(){
		boolean res = false;
		if(found_eFT && fkps.size() >= 2 && ckms.size() >= 1 ){
			double avg_vz = (double)(fkps.get(0).vz() + fkps.get(1).vz() + ckms.get(0).vz())/3;
			Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Vftel_corrected.vector());Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(fkps.get(1).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			efkpfkpckm_MM_efkpfkp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Vftel_corrected.vector());Vmissekpkpkm.sub(fkps.get(0).vector());Vmissekpkpkm.sub(fkps.get(1).vector());Vmissekpkpkm.sub(ckms.get(0).vector());
			efkpfkpckm_MM_efkpfkpckm = (float)Vmissekpkpkm.mass();

			if (Math.abs(efkpfkpckm_MM_efkpfkpckm - 1.115683) <= 0.052){
				efkpfkpckm_found_lambda = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedlambda = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.115683f*1.115683f));
				LorentzVector VxifromLambda = new LorentzVector(0, 0, 0, 0);
				VxifromLambda.add(Vconstrainedlambda);VxifromLambda.add(ckms.get(0).vector());
				efkpfkpckm_IM_kmlambda = (float)VxifromLambda.mass();
				
			}

			res = true;
		}
		return res;
	}

	public boolean select_efkpckpckm(){
		boolean res = false;
		if(found_eFT && fkps.size() >= 1&& ckps.size()>=1 && ckms.size() >= 1 ){
			double avg_vz = (double)(fkps.get(0).vz() + ckps.get(0).vz() + ckms.get(0).vz())/3;
			Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Vftel_corrected.vector());Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			efkpckpckm_MM_efkpckp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Vftel_corrected.vector());Vmissekpkpkm.sub(fkps.get(0).vector());Vmissekpkpkm.sub(ckps.get(0).vector());Vmissekpkpkm.sub(ckms.get(0).vector());
			efkpckpckm_MM_efkpckpckm = (float)Vmissekpkpkm.mass();

			if (Math.abs(efkpckpckm_MM_efkpckpckm - 1.115683) <= 0.052){
				efkpckpckm_found_lambda = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedlambda = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.115683f*1.115683f));
				LorentzVector VxifromLambda = new LorentzVector(0, 0, 0, 0);
				VxifromLambda.add(Vconstrainedlambda);VxifromLambda.add(ckms.get(0).vector());
				efkpckpckm_IM_kmlambda = (float)VxifromLambda.mass();
				
			}

			res = true;
		}
		return res;
	}

	public boolean select_eckpckpckm(){
		boolean res = false;
		if(found_eFT && ckps.size() >= 2 && ckms.size() >= 1 ){
			double avg_vz = (double)(ckps.get(0).vz() + ckps.get(1).vz() + ckms.get(0).vz())/3;
			Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Vftel_corrected.vector());Vmissekpkp.sub(ckps.get(0).vector());Vmissekpkp.sub(ckps.get(1).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			eckpckpckm_MM_eckpckp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Vftel_corrected.vector());Vmissekpkpkm.sub(ckps.get(0).vector());Vmissekpkpkm.sub(ckps.get(1).vector());Vmissekpkpkm.sub(ckms.get(0).vector());
			eckpckpckm_MM_eckpckpckm = (float)Vmissekpkpkm.mass();

			if (Math.abs(eckpckpckm_MM_eckpckpckm - 1.115683) <= 0.052){
				eckpckpckm_found_lambda = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedlambda = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.115683f*1.115683f));
				LorentzVector VxifromLambda = new LorentzVector(0, 0, 0, 0);
				VxifromLambda.add(Vconstrainedlambda);VxifromLambda.add(ckms.get(0).vector());
				eckpckpckm_IM_kmlambda = (float)VxifromLambda.mass();
				
			}

			res = true;
		}
		return res;
	}

	public boolean select_efkpckp(){
		boolean res = false;
		if((found_eFT && fkps.size() == 1 && ckps.size()==1)) {
			res = true;
			float fkpckp_deltap = fkps.get(0).p()-ckps.get(0).p();
			float fkpckp_deltatheta =  Math.toDegrees(fkps.get(0).theta())-Math.toDegrees(ckps.get(0).theta());
			float fkpckp_deltaphi = Math.toDegrees(fkps.get(0).phi())-Math.toDegrees(ckps.get(0).phi());
			//if((fkpckp_deltap > -0.5 || fkpckp_deltap < 0.1) && (fkpckp_deltatheta > -8.0) && (fkpckp_deltaphi > -25.0 || fkpckp_deltaphi < 10.0)) res = true;
		}
	//	/*
		else if(found_eFT && kps.size() == 1 && fprots.size() == 1 && cprots.size()==1){
			res = true;
			float fkpckp_deltap = fprots.get(0).p()-cprots.get(0).p();
			float fkpckp_deltatheta =  Math.toDegrees(fprots.get(0).theta())-Math.toDegrees(cprots.get(0).theta());
			float fkpckp_deltaphi = Math.toDegrees(fprots.get(0).phi())-Math.toDegrees(cprots.get(0).phi());
			//if((fkpckp_deltap > -0.5 || fkpckp_deltap < 0.1) && (fkpckp_deltatheta > -8.0) && (fkpckp_deltaphi > -25.0 || fkpckp_deltaphi < 10.0)) res = true;
		}
		  /*
		else if(found_eFT && kps.size() == 1 && fpips.size() == 1 && cpips.size()==1) {
			res = true;
			float fkpckp_deltap = fpips.get(0).p()-cpips.get(0).p();
			float fkpckp_deltatheta =  Math.toDegrees(fpips.get(0).theta())-Math.toDegrees(cpips.get(0).theta());
			float fkpckp_deltaphi = Math.toDegrees(fpips.get(0).phi())-Math.toDegrees(cpips.get(0).phi());
			//if((fkpckp_deltap > -0.5 || fkpckp_deltap < 0.1) && (fkpckp_deltatheta > -8.0) && (fkpckp_deltaphi > -25.0 || fkpckp_deltaphi < 10.0)) res = true;
		}
		//|| (found_eFT && fpips.size() >= 1 && cpips.size()>=1) || (found_eFT && fprots.size() >= 1 && cprots.size()>=1)){
			//res = true;
		//} */

		return res;
	}




	public boolean select_ekpkp_req_km(){
		boolean res = false;
		//ftels = new ArrayList<Particle>();
		//if( found_eFT && kps.size() >= 2 && kms.size() >= 1 && mcels.size() ==1 && mckps.size() == 2 && mckms.size() == 1) {
		if( found_eFT && kps.size() == 2 && kms.size() == 1 ) {
			// labeling generated kps with momentum
			/*if (mckps.get(0).p() > mckps.get(1).p()){
				mckp1 = mckps.get(0);
				mckp2 = mckps.get(1);
			} else {
				mckp1 = mckps.get(1);
				mckp2 = mckps.get(0);
			}*/

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
			//LorentzVector Ve_corrected = electron_corn(Ve, e_ftCal_hitPosition, avg_vz);
			Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			//Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, (float)mcels.get(0).vz());
			// to use generated phi for the detected electron to see if one can achieve improvements with phi correction
			//Vftel_corrected.vector().vect().setMagThetaPhi(Vftel_corrected.vector().p(), Vftel_corrected.vector().theta(), (float)mcels.get(0).vector().phi());
			//Vftel_corrected.vector().vect().setMagThetaPhi(Vftel_corrected.vector().p(), (float)mcels.get(0).vector().theta(), Vftel_corrected.vector().phi());
			//Vftel_corrected.vector().vect().setMagThetaPhi(Vftel_corrected.vector().p(), (float)mcels.get(0).vector().theta(), (float)mcels.get(0).vector().phi());
			//LorentzVector Ve_corrected = electron_corn(Ve, e_ftCal_hitPosition, mcels.get(0).vz());
			LorentzVector Vmissekpkp_noCorr = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp_noCorr.add(VT);Vmissekpkp_noCorr.add(VB);Vmissekpkp_noCorr.sub(Ve);Vmissekpkp_noCorr.sub(reckp1.vector());Vmissekpkp_noCorr.sub(reckp2.vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			ekpkpkm_MM_ekpkp_nocorr = (float)Vmissekpkp_noCorr.mass();

			LorentzVector Vmissekpkpkm_noCorr = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm_noCorr.add(VT);Vmissekpkpkm_noCorr.add(VB);Vmissekpkpkm_noCorr.sub(Ve);Vmissekpkpkm_noCorr.sub(reckp1.vector());Vmissekpkpkm_noCorr.sub(reckp2.vector());Vmissekpkpkm_noCorr.sub(kms.get(0).vector());
			ekpkpkm_MM_ekpkpkm_nocorr = (float)Vmissekpkpkm_noCorr.mass();


			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Vftel_corrected.vector());Vmissekpkp.sub(reckp1.vector());Vmissekpkp.sub(reckp2.vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			ekpkpkm_MM_ekpkp = (float)Vmissekpkp.mass();



			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Vftel_corrected.vector());Vmissekpkpkm.sub(reckp1.vector());Vmissekpkpkm.sub(reckp2.vector());Vmissekpkpkm.sub(kms.get(0).vector());
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
	public Particle electron_energy_correction(Particle el){
		Double e_new, px_el, py_el, pz_el;
		
		e_new = el.e()+0.1574-0.01178*el.e()-0.007345*el.e()*el.e()+0.002909*el.e()*el.e()*el.e()-0.0002243*el.e()*el.e()*el.e()*el.e();
		px_el = e_new*(el.px()/el.p());
		py_el = e_new*(el.py()/el.p());
		pz_el = e_new*(el.pz()/el.p());
		Particle el_new = new Particle(11, 0, 0, 0, 0, 0, 0);
		el_new.setVector(11, px_el, py_el, pz_el, el.vx(), el.vy(), el.vz());
		return el_new;
	}




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
		efkpfkpfkm_found_lambda = false;
		efkpckpfkm_found_lambda = false;
		eckpckpfkm_found_lambda = false;
		efkpfkpckm_found_lambda = false;
		efkpckpckm_found_lambda = false;
		eckpckpckm_found_lambda = false;
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
		reckp1 = null;
		reckp2 = null;
		/*
		pips.clear();
		kps.clear();
		kms.clear();
		pims.clear();
		prots.clear();
		ckps.clear();
		fkps.clear();
		fpips.clear();
		fkms.clear();
		fpims.clear();
		fprots.clear();
		ckps.clear();
		cpips.clear();
		ckms.clear();
		cpims.clear();
		cprots.clear();
	//	*/
		/*
		pips = null;
		kps = null;
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
	//	*/
			
	}

	public void processEvent(DataEvent event) {
		resetCounters();

		//Timothy's analysis fitter
        //GenericKinematicFitter research_fitter = new analysis_fitter(Eb);

				/*
				double eb = 10.604;
				analysis_fitter af = new analysis_fitter(eb);
				if (af.bank_test(event)){
					System.out.println("PassedEvent = " + event);
				}
				//*/



		if (event.hasBank("RECFT::Event"))
			fillRecBank(event.getBank("RECFT::Event"));
		if (event.hasBank("REC::Particle") && event.hasBank("RECFT::Particle") && event.hasBank("REC::ForwardTagger")) e_ft_part_ind = makeFTElectron(event.getBank("REC::Particle"), event.getBank("RECFT::Particle"), event.getBank("REC::ForwardTagger"));
		if (event.hasBank("MC::Particle") == true ) fillMCPartBank(event.getBank("MC::Particle"), event.getBank("RUN::config"));
		if(e_ft_part_ind > -1 && found_eFT) {		
			if (event.hasBank("REC::Particle") && event.hasBank("RECFT::Particle") && event.hasBank("REC::Scintillator")) makeOthers(event.getBank("REC::Particle"), event.getBank("RECFT::Particle"), event.getBank("REC::Scintillator"));
			
			//if (event.hasBank("REC::Scintillator")) fillFTOF(event.getBank("REC::Particle"), event.getBank("REC::Scintillator"));
			//System.out.println("found pip :: " + pips.size() + " pips tracks in FD & " + cpips.size() + " pip tracks in CD.");
			FillHists();
			
		} // e_ft_part_ind > -1
		
		
	} //processEvent

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
    		f1.setParameter(2, 0.1);
    		DataFitter.fit(f1, hslice.get(i), "Q"); //No options uses error for sigma 
    		if(amp>10) meanX.addPoint(x, f1.getParameter(1), ex, f1.parameter(1).error());
		}
	}

	public void plotGraph(){
		// for tracks reconstructed twice in FD and CD
		H2F h1 = dg_fkpckp_pthphi.getH2F("hi-fkpckp-deltap-p");
		//H2F h1 = hi_fkpckp_deltap_p;
		F1D f1 = new F1D("f1","[amp]*gaus(x,[mean],[sigma])", -1.0, 1.0);
		GraphErrors gr1 = new GraphErrors();
		sliceFit(h1, f1, gr1);
		myCanvas.getCanvas("fkpckp").cd(7);
		myCanvas.getCanvas("fkpckp").draw(gr1);

		H2F h2 = dg_fkpckp_pthphi.getH2F("hi-fkpckp-deltap-theta");
		//H2F h2 = hi_fkpckp_deltap_theta;
		F1D f2 = new F1D("f2","[amp]*gaus(x,[mean],[sigma])", -1.0, 1.0);
		GraphErrors gr2 = new GraphErrors();
		sliceFit(h2, f2, gr2);
		myCanvas.getCanvas("fkpckp").cd(11);
		myCanvas.getCanvas("fkpckp").draw(gr2);

		H2F h3 = dg_fkpckp_pthphi.getH2F("hi-fkpckp-deltap-phi");
		//H2F h3 = hi_fkpckp_deltap_phi
		F1D f3 = new F1D("f3","[amp]*gaus(x,[mean],[sigma])", -1.0, 1.0);
		GraphErrors gr3 = new GraphErrors();
		sliceFit(h3, f3, gr3);
		myCanvas.getCanvas("fkpckp").cd(15);
		myCanvas.getCanvas("fkpckp").draw(gr3);

	}

	
	public void analyze() {

		/*
		fitxi(dg_rec_xi.getH1F("hi_ekpkpkm_MM_ekpkp_nocorr"), dg_rec_xi.getF1D("fn_xi_no_corr"));
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
    	
		*/

    	//fitlambda(dg_rec_xi.getH1F("hi_ekpkpkm_MM_ekpkpkm"), dg_rec_xi.getF1D("f1_lambda"));
    	//fitlambda(dg_rec_xi.getH1F("hi_ekpkpkm_MM_ekpkpkm"), dg_rec_xi.getF1D("f1_sigma"));
    	//fitvz(dg_rec_kp1.getH1F("hi_rec_kp1_dvz"), dg_rec_kp1.getF1D("f_gaus"));
    	//fitvz(dg_rec_kp2.getH1F("hi_rec_kp2_dvz"), dg_rec_kp2.getF1D("f_gaus"));
    	//fitvz(dg_rec_km.getH1F("hi_rec_km_dvz"), dg_rec_km.getF1D("f_gaus"));
    //	fitxi(dg_mm.getH1F("hi_"))


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
//			hi_FT_e_t_f.fill(e_phi, e_the);
//			hi_FT_e_p_f.fill(e_phi, e_mom);
//			hi_FT_e_p_the.fill(e_the, e_mom);
		/*
			hi_FT_W_Q2.fill(e_W, e_Q2);
			hi_FT_W.fill(e_W);
			hi_FT_Q2.fill(e_Q2);	
			hi_virphoton.fill(e_virphoton);	
			*/
		}

		// fill particle counters with one electron detected in FT

		
		if(select_ekpkp_req_pim()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			//hi_ekpkp_MM_req_pim.fill(ekpkp_MM_req_pim);
		}
		if(select_efkpckp_req_pim()) {
			/*
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			hi_efkpckp_MM_req_pim.fill(efkpckp_MM_req_pim);
			*/
		}
		
		if(select_efkpckp_req_cpim()) {
			/*
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			hi_efkpckp_MM_req_cpim.fill(efkpckp_MM_req_cpim);
			*/
		}
		
		if(select_efkpckp_req_km()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			hi_efkpckp_MM_req_km.fill(efkpckp_MM_req_km);
		}
		
		if(select_efkpckp_req_ckm()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			hi_efkpckp_MM_req_ckm.fill(efkpckp_MM_req_ckm);
		}

		if(select_ekpkpkmprot()){
			hi_ekpkpkmprot_MM2.fill(ekpkpkmprot_MM2);
		}

		// ekpkpkm required 
		if(select_efkpfkpfkm() && !select_efkpckpfkm() && !select_efkpfkpckm()){

			hi_FT_W_Q2.fill(e_W, e_Q2);
			hi_FT_W.fill(e_W);
			hi_FT_Q2.fill(e_Q2);	
			hi_virphoton.fill(e_virphoton);

			hi_efkpfkpfkm_MM_efkpfkp.fill(efkpfkpfkm_MM_efkpfkp);
			hi_efkpfkpfkm_MM_efkpfkpfkm.fill(efkpfkpfkm_MM_efkpfkpfkm);
			hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm.fill(efkpfkpfkm_MM_efkpfkpfkm, efkpfkpfkm_MM_efkpfkp);

			if(efkpfkpfkm_found_lambda) {
				hi_efkpfkpfkm_MM_efkpfkp_lam_evnt.fill(efkpfkpfkm_MM_efkpfkp);
				hi_efkpfkpfkm_IM_kmlambda.fill(efkpfkpfkm_IM_kmlambda);
			}

			hi_ekpkpkm_MM_ekpkp.fill(efkpfkpfkm_MM_efkpfkp);
			hi_ekpkpkm_MM_ekpkpkm.fill(efkpfkpfkm_MM_efkpfkpfkm);


			// labeling reconstructed kps with momentum
			Particle Vslowkp;
			Particle Vfastkp;

			if (fkps.get(0).p() > fkps.get(1).p()){
				Vfastkp = new Particle(fkps.get(0));
				Vslowkp = new Particle(fkps.get(1));
			} else {
				Vfastkp = new Particle(fkps.get(1));
				Vslowkp = new Particle(fkps.get(0));
			}

			Particle p_fte;

			if (runType == "data"){

				// Applying simple geometrical correction for energy corrected FT electron. Geometrical correction derived based on common vertex position of 3 kaons look:: electron_corn(Particle(), Vector3D(), double)
				double avg_vz = (double)(fkps.get(0).vz() + fkps.get(1).vz() + fkms.get(0).vz())/3;
				p_fte = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz); // geometrical correction on top of ft energy correction
				//p_fte = ftels.get(0); //FT energy correction only
				//p_fte = Velectron; //no corn at all
			} else if (runType == "mc"){
				p_fte = Velectron; // for the simulation we don't need to apply data driven FT electron energy correction and geometrical correction
			} else {
                println("You provided wrong argument for runType. Please specify mc or data")
            }


			// Applying simple geometrical correction for energy corrected FT electron. Geometrical correction derived based on common vertex position of 3 kaons look:: electron_corn(Particle(), Vector3D(), double)
			//double avg_vz = (double)(fkps.get(0).vz() + fkps.get(1).vz() + fkms.get(0).vz())/3;
			//p_fte = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);

			//p_fte = ftels.get(0);
			//p_fte = Velectron; // for the simulation we don't need to apply data driven FT electron energy correction and geometrical correction

			// set up boost to gamma*-nucleon center of mass frame
			LorentzVector lv_q = new LorentzVector(VB); lv_q.sub(p_fte.vector());
			//LorentzVector lv_q = new LorentzVector(VB); lv_q.sub(p_fte.vector()); // For Simulation Try without applying correction to electron
        	LorentzVector gN = new LorentzVector(lv_q);
			gN.add(VT);
			Vector3 gNBoost = gN.boostVector();
			gNBoost.negative();

			// boost to gamma*-nucleon center of mass frame
			LorentzVector lv_target_gN = new LorentzVector(VT); lv_target_gN.boost(gNBoost);
			LorentzVector lv_e_gN = new LorentzVector(p_fte.vector()); lv_e_gN.boost(gNBoost);
			LorentzVector lv_fastkp_gN = new LorentzVector(Vfastkp.vector()); lv_fastkp_gN.boost(gNBoost);
			LorentzVector lv_slowkp_gN = new LorentzVector(Vslowkp.vector()); lv_slowkp_gN.boost(gNBoost);
			LorentzVector lv_km_gN = new LorentzVector(fkms.get(0).vector()); lv_km_gN.boost(gNBoost);

			//for mass calculation
			LorentzVector lv_fkpskp = new LorentzVector();lv_fkpskp.add(Vfastkp.vector());lv_fkpskp.add(Vslowkp.vector());
			double m_fkpskp = lv_fkpskp.mass();
			LorentzVector lv_fkpskpkm = new LorentzVector();lv_fkpskpkm.add(lv_fkpskp);lv_fkpskpkm.add(fkms.get(0).vector());
			double m_fkpskpkm = lv_fkpskpkm.mass();
			LorentzVector lv_fkpkm = new LorentzVector();lv_fkpkm.add(Vfastkp.vector());lv_fkpkm.add(fkms.get(0).vector());
			double m_fkpkm = lv_fkpkm.mass();
			LorentzVector lv_skpkm = new LorentzVector();lv_skpkm.add(Vslowkp.vector());lv_fkpkm.add(fkms.get(0).vector());
			double m_skpkm = lv_skpkm.mass();
			//missing mass calculation
			LorentzVector lv_mm_e = new LorentzVector();lv_mm_e.add(lv_q);lv_mm_e.add(VT);
			double mm_e = lv_mm_e.mass();
			LorentzVector lv_mm_fkp = new LorentzVector();lv_mm_fkp.add(lv_mm_e);lv_mm_fkp.sub(Vfastkp.vector());
			double mm_fkp = lv_mm_fkp.mass();
			LorentzVector lv_mm_skp = new LorentzVector();lv_mm_skp.add(lv_mm_e);lv_mm_skp.sub(Vslowkp.vector());
			double mm_skp = lv_mm_skp.mass();
			LorentzVector lv_mm_km = new LorentzVector();lv_mm_km.add(lv_mm_e);lv_mm_km.sub(fkms.get(0).vector());
			double mm_km =  lv_mm_km.mass();
			LorentzVector lv_mm_fkpskp = new LorentzVector();lv_mm_fkpskp.add(lv_mm_fkp);lv_mm_fkpskp.sub(Vslowkp.vector());
			double mm_fkpskp = lv_mm_fkpskp.mass();
			LorentzVector lv_mm_fkpkm = new LorentzVector();lv_mm_fkpkm.add(lv_mm_fkp);lv_mm_fkpkm.sub(fkms.get(0).vector());
			double mm_fkpkm = lv_mm_fkpkm.mass();
			LorentzVector lv_mm_skpkm = new LorentzVector();lv_mm_skpkm.add(lv_mm_skp);lv_mm_skpkm.sub(fkms.get(0).vector());
			double mm_skpkm = lv_mm_skpkm.mass();
			LorentzVector lv_mm_fkpskpkm = new LorentzVector();lv_mm_fkpskpkm.add(lv_mm_fkpskp);lv_mm_fkpskpkm.sub(fkms.get(0).vector());
			double mm_fkpskpkm = lv_mm_fkpskpkm.mass();

			//println("MM: "+ mm_fkpskpkm);

			/*
			//Store three kaon event information in the arraylist so that we can implement event mixing technique later on
            //public static ArrayList<LorentzVector> lv_my_ftels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
            lv_my_ftels.add(p_fte.vector());
            lv_my_fast_fkps.add(Vfastkp.vector());
            lv_my_slow_fkps.add(Vslowkp.vector());
            lv_my_fkms.add(fkms.get(0).vector());
            //*/

			///*
			//public static ArrayList<LorentzVector> lv_my_xi_ftels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms;
			//Xi data sample from MM(eK+K+K-) cut to select lambda/sigma region from the fit and store information so tha we can implement event mixing technique later on to estimate bg shape
			if (1.00 <= mm_fkpskpkm && mm_fkpskpkm <= 1.29){ //f18out:: 1.00 <= mm_fkpskpkm && mm_fkpskpkm <= 1.29 // f18in:: 1.02 <= mm_fkpskpkm && mm_fkpskpkm <= 1.28// s19:: 1.01 <= mm_fkpskpkm && mm_fkpskpkm <= 1.27
				lv_my_xi_ftels.add(p_fte.vector());
            	lv_my_xi_fast_fkps.add(Vfastkp.vector());
            	lv_my_xi_slow_fkps.add(Vslowkp.vector());
            	lv_my_xi_fkms.add(fkms.get(0).vector());
			}
			//*/

			// considering fast kp as pip
			LorentzVector lv_fkpaspip = new LorentzVector();lv_fkpaspip.setPxPyPzM(Vfastkp.px(), Vfastkp.py(), Vfastkp.pz(), PDGDatabase.getParticleById(211).mass());
			LorentzVector lv_skpaspip = new LorentzVector();lv_skpaspip.setPxPyPzM(Vslowkp.px(), Vslowkp.py(), Vslowkp.pz(), PDGDatabase.getParticleById(211).mass());
			LorentzVector lv_kmaspim = new LorentzVector();lv_kmaspim.setPxPyPzM(fkms.get(0).px(), fkms.get(0).py(), fkms.get(0).pz(), PDGDatabase.getParticleById(211).mass());

			//LorentzVector lv_mm_fkpaspipskpkm = new LorentzVector();lv_mm_fkpaspipskpkm.add(VT);lv_mm_fkpaspipskpkm.add(VB);lv_mm_fkpaspipskpkm.sub(Ve);lv_mm_fkpaspipskpkm.sub(lv_fkpaspip);lv_mm_fkpaspipskpkm.sub(Vslowkp.vector());lv_mm_fkpaspipskpkm.sub(fkms.get(0).vector());
			LorentzVector lv_mm_fkpaspipskpkm = new LorentzVector();lv_mm_fkpaspipskpkm.add(lv_mm_skpkm);lv_mm_fkpaspipskpkm.sub(lv_fkpaspip);
			//hi_efkpfkpfkm_MM_efastkpaspipslowkpfkm.fill(lv_mm_fkpaspipskpkm.mass(), efkpfkpfkm_MM_efkpfkpfkm);
			double mm_fkpaspipskpkm = lv_mm_fkpaspipskpkm.mass();

			//// considering slow kp as pip
			LorentzVector lv_mm_fkpskpaspipkm = new LorentzVector();lv_mm_fkpskpaspipkm.add(lv_mm_fkpkm);lv_mm_fkpskpaspipkm.sub(lv_skpaspip);
			//hi_efkpfkpfkm_MM_eslowkpaspipfastkpfkm.fill(lv_mm_fkpskpaspipkm.mass(), efkpfkpfkm_MM_efkpfkpfkm);
			double mm_fkpskpaspipkm = lv_mm_fkpskpaspipkm.mass();

			LorentzVector lv_mm_fkpskpkmaspim = new LorentzVector();lv_mm_fkpskpkmaspim.add(lv_mm_fkpskp);lv_mm_fkpskpkmaspim.sub(lv_kmaspim);
			double mm_fkpskpkmaspim = lv_mm_fkpskpkmaspim.mass();
			//LorentzVector lv_mm_fkpaspipskpaspipkm = new LorentzVector();lv_mm_fkpaspipskpaspipkm.add(VT);lv_mm_fkpaspipskpaspipkm.add(VB);lv_mm_fkpaspipskpaspipkm.sub(Ve);lv_mm_fkpaspipskpaspipkm.sub(lv_fkpaspip);lv_mm_fkpaspipskpaspipkm.sub(lv_skpaspip);lv_mm_fkpaspipskpaspipkm.sub(lv_kmaspim);
			LorentzVector lv_mm_fkpaspipskpaspipkm = new LorentzVector();lv_mm_fkpaspipskpaspipkm.add(lv_mm_km);lv_mm_fkpaspipskpaspipkm.sub(lv_fkpaspip);lv_mm_fkpaspipskpaspipkm.sub(lv_skpaspip);
			double mm_fkpaspipskpaspipkm = lv_mm_fkpaspipskpaspipkm.mass();

			LorentzVector lv_mm_fkpaspipskpkmaspim = new LorentzVector();lv_mm_fkpaspipskpkmaspim.add(lv_mm_skp);lv_mm_fkpaspipskpkmaspim.sub(lv_fkpaspip);lv_mm_fkpaspipskpkmaspim.sub(lv_kmaspim);
			double mm_fkpaspipskpkmaspim = lv_mm_fkpaspipskpkmaspim.mass();

			LorentzVector lv_mm_fkpskpaspipkmaspim = new LorentzVector();lv_mm_fkpskpaspipkmaspim.add(lv_mm_fkp);lv_mm_fkpskpaspipkmaspim.sub(lv_skpaspip);lv_mm_fkpskpaspipkmaspim.sub(lv_kmaspim);
			double mm_fkpskpaspipkmaspim = lv_mm_fkpskpaspipkmaspim.mass();

			LorentzVector lv_mm_fkpaspipskpaspipkmaspim = new LorentzVector();lv_mm_fkpaspipskpaspipkmaspim.add(lv_mm_e);lv_mm_fkpaspipskpaspipkmaspim.sub(lv_fkpaspip);lv_mm_fkpaspipskpaspipkmaspim.sub(lv_skpaspip);lv_mm_fkpaspipskpaspipkmaspim.sub(lv_kmaspim);
			double mm_fkpaspipskpaspipkmaspim = lv_mm_fkpaspipskpaspipkmaspim.mass();

			//electron kinematics
			double q2 = -lv_q.mass2();
			double nu = VB.e()-p_fte.vector().e();
			//double nu = VB.e()-p_fte.vector().e(); //for simulation try without applying correction to electron
			double x  = q2 / (2 * PDGDatabase.getParticleById(2212).mass() * nu);
			double W  = Math.pow(Math.pow(PDGDatabase.getParticleById(2212).mass(),2)+2*PDGDatabase.getParticleById(2212).mass()*nu - q2, 0.5);

			// lab e kinematics
			double beam_e = VB.e();
			///*
			double e_px = p_fte.vector().px();
			double e_py = p_fte.vector().py();
			double e_pz = p_fte.vector().pz();
			double e_p = p_fte.vector().p();
			double e_e = p_fte.vector().e();
			double e_vx = p_fte.vx();
			double e_vy = p_fte.vy();
			double e_vz = p_fte.vz();
			//*/

			// For simulation try without correction change p_fte to ftels.get(0) {It has FT electron energy correction only} or Velectron it has no corrn at all
		/*
			double e_px = p_fte.vector().px();
			double e_py = p_fte.vector().py();
			double e_pz = p_fte.vector().pz();
			double e_p = p_fte.vector().p();
			double e_e = p_fte.vector().e();
			double e_vx = p_fte.vx();
			double e_vy = p_fte.vy();
			double e_vz = p_fte.vz();
			*/

			double e_theta = Math.toDegrees(p_fte.vector().theta());
			double e_phi = Math.toDegrees(p_fte.vector().phi());
			// lab fkp kinematics
			double fkp_px = Vfastkp.px();
			double fkp_py = Vfastkp.py();
			double fkp_pz = Vfastkp.pz();
			double fkp_p = Vfastkp.p();
			double fkp_e = Vfastkp.e();
			double fkp_vx = Vfastkp.vx();
			double fkp_vy = Vfastkp.vy();
			double fkp_vz = Vfastkp.vz();
			double fkp_theta = Math.toDegrees(Vfastkp.theta());
			double fkp_phi = Math.toDegrees(Vfastkp.phi());
			// lab skp kinematics
			double skp_px = Vslowkp.px();
			double skp_py = Vslowkp.py();
			double skp_pz = Vslowkp.pz();
			double skp_p = Vslowkp.p();
			double skp_e = Vslowkp.e();
			double skp_vx = Vslowkp.vx();
			double skp_vy = Vslowkp.vy();
			double skp_vz = Vslowkp.vz();
			double skp_theta = Math.toDegrees(Vslowkp.theta());
			double skp_phi = Math.toDegrees(Vslowkp.phi());
			// lab km kinematics
			double km_px = fkms.get(0).px();
			double km_py = fkms.get(0).py();
			double km_pz = fkms.get(0).pz();
			double km_p = fkms.get(0).p();
			double km_e = fkms.get(0).e();
			double km_vx = fkms.get(0).vx();
			double km_vy = fkms.get(0).vy();
			double km_vz = fkms.get(0).vz();
			double km_theta = Math.toDegrees(fkms.get(0).theta());
			double km_phi = Math.toDegrees(fkms.get(0).phi());

			// append event to next line of the text file
			//file.append(runnum+" "+evnum+" "+helicity+" ");
			outFile.append(runnum+" "+evnum+" "+beam_e+" "+q2+" "+nu+" "+x+" "+W+" ");
			outFile.append(e_px+" "+e_py+" "+e_pz+" "+e_p+" "+e_e+" "+e_vx+" "+e_vy+" "+e_vz+" "+e_theta+" "+e_phi+" ");
			outFile.append(fkp_px+" "+fkp_py+" "+fkp_pz+" "+fkp_p+" "+fkp_e+" "+fkp_vx+" "+fkp_vy+" "+fkp_vz+" "+fkp_theta+" "+fkp_phi+" ");
			outFile.append(skp_px+" "+skp_py+" "+skp_pz+" "+skp_p+" "+skp_e+" "+skp_vx+" "+skp_vy+" "+skp_vz+" "+skp_theta+" "+skp_phi+" ");
			outFile.append(km_px+" "+km_py+" "+km_pz+" "+km_p+" "+km_e+" "+km_vx+" "+km_vy+" "+km_vz+" "+km_theta+" "+km_phi+" ");
			outFile.append(m_fkpskp+" "+m_fkpkm+" "+m_skpkm+" "+m_fkpskpkm+" ");
			outFile.append(mm_e+" "+mm_fkp+" "+mm_skp+" "+mm_km+" "+mm_fkpskp+" "+mm_fkpkm+" "+mm_skpkm+" "+mm_fkpskpkm+" ");
			outFile.append(mm_fkpaspipskpkm+" "+mm_fkpskpaspipkm+" "+mm_fkpskpkmaspim+" "+mm_fkpaspipskpaspipkm+" "+mm_fkpaspipskpkmaspim+" "+mm_fkpskpaspipkmaspim+" "+mm_fkpaspipskpaspipkmaspim+"\n");
			//outFile.append("\n");
			//println(); println();
			//print("1:runnum, 2:evnum, 3:beam_e, 4:q2, 5:nu, 6:x, 7:W, 8:e_px, 9:e_py, 10:e_pz, 11:e_p, 12:e_e, 13:e_vx, 14:e_vy, 15:e_vz, 16:e_theta, 17:e_phi,");
			//print("14:fkp_px, 15:fkp_py, 16:fkp_pz, 17:fkp_p, 18:fkp_e, 19:fkp_vx, 20:fkp_vy, 21:fkp_vz, 22:fkp_theta, 23:fkp_phi,");
			//print("24:skp_px, 25:skp_py, 26:skp_pz, 27:skp_p, 28:skp_e, 29:skp_vx, 30:skp_vy, 31:skp_vz, 32:skp_theta, 33:skp_phi,");
			//print("34:km_px, 35:km_py, 36:km_pz, 37:km_p, 38:km_e, 39:km_vx, 40:km_vy, 41:km_vz, 42:km_theta, 43:km_phi,");






		}
		if(select_efkpckpfkm() && !select_efkpfkpfkm() && !select_efkpfkpckm()){
			float fkpckp_deltap = fkps.get(0).p()-ckps.get(0).p();
			float fkpckp_deltatheta =  Math.toDegrees(fkps.get(0).theta())-Math.toDegrees(ckps.get(0).theta());
			float fkpckp_deltaphi = Math.toDegrees(fkps.get(0).phi())-Math.toDegrees(ckps.get(0).phi());
			if ((fkpckp_deltap < -0.5 || fkpckp_deltap > 0.1) && (fkpckp_deltatheta < -8.0) && (fkpckp_deltaphi < -25.0 || fkpckp_deltaphi > 10.0) ){
				hi_efkpckpfkm_MM_efkpckpfkm_nocorn.fill(efkpckpfkm_MM_efkpckpfkm_nocorn);
				hi_efkpckpfkm_MM_efkpckp.fill(efkpckpfkm_MM_efkpckp);
				hi_efkpckpfkm_MM_efkpckpfkm.fill(efkpckpfkm_MM_efkpckpfkm);
				hi_efkpckpfkm_MM_efkpckp_MM_efkpckpfkm.fill(efkpckpfkm_MM_efkpckpfkm, efkpckpfkm_MM_efkpckp);

				if(efkpckpfkm_found_lambda){
					hi_efkpckpfkm_MM_efkpckp_lam_evnt.fill(efkpckpfkm_MM_efkpckp);
					hi_efkpckpfkm_IM_kmlambda.fill(efkpckpfkm_IM_kmlambda);
				}
			}
		}
		if(select_eckpckpfkm()){
			hi_eckpckpfkm_MM_eckpckp.fill(eckpckpfkm_MM_eckpckp);
			hi_eckpckpfkm_MM_eckpckpfkm.fill(eckpckpfkm_MM_eckpckpfkm);
			hi_eckpckpfkm_MM_eckpckp_MM_eckpckpfkm.fill(eckpckpfkm_MM_eckpckpfkm, eckpckpfkm_MM_eckpckp);
			if(eckpckpfkm_found_lambda){
				hi_eckpckpfkm_MM_eckpckp_lam_evnt.fill(eckpckpfkm_MM_eckpckp);
				hi_eckpckpfkm_IM_kmlambda.fill(eckpckpfkm_IM_kmlambda);
			}
		}
		if(select_efkpfkpckm()&& !select_efkpckpfkm() && !select_efkpfkpfkm()){
			hi_efkpfkpckm_MM_efkpfkp.fill(efkpfkpckm_MM_efkpfkp);
			hi_efkpfkpckm_MM_efkpfkpckm.fill(efkpfkpckm_MM_efkpfkpckm);
			hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm.fill(efkpfkpckm_MM_efkpfkpckm, efkpfkpckm_MM_efkpfkp);
			if(efkpfkpckm_found_lambda){
				hi_efkpfkpckm_MM_efkpfkp_lam_evnt.fill(efkpfkpckm_MM_efkpfkp);
				hi_efkpfkpckm_IM_kmlambda.fill(efkpfkpckm_IM_kmlambda);
			}
		}
		if(select_efkpckpckm()){
			float fkpckp_deltap = fkps.get(0).p()-ckps.get(0).p();
			float fkpckp_deltatheta =  Math.toDegrees(fkps.get(0).theta())-Math.toDegrees(ckps.get(0).theta());
			float fkpckp_deltaphi = Math.toDegrees(fkps.get(0).phi())-Math.toDegrees(ckps.get(0).phi());
			//hi_fkpckp_deltap.fill(fkpckp_deltap);
			//hi_fkpckp_deltatheta.fill(fkpckp_deltatheta);
			//hi_fkpckp_deltaphi.fill(fkpckp_deltaphi);
			if ((fkpckp_deltap < -0.5 || fkpckp_deltap > 0.1) && (fkpckp_deltatheta < -8.0) && (fkpckp_deltaphi < -25.0 || fkpckp_deltaphi > 10.0) ){
				hi_efkpckpckm_MM_efkpckp.fill(efkpckpckm_MM_efkpckp);
				hi_efkpckpckm_MM_efkpckpckm.fill(efkpckpckm_MM_efkpckpckm);
				hi_efkpckpckm_MM_efkpckp_MM_efkpckpckm.fill(efkpckpckm_MM_efkpckpckm, efkpckpckm_MM_efkpckp);
				if(efkpckpckm_found_lambda){
					hi_efkpckpckm_MM_efkpckp_lam_evnt.fill(efkpckpckm_MM_efkpckp);
					hi_efkpckpckm_IM_kmlambda.fill(efkpckpckm_IM_kmlambda);
				}
			}
			
		}

		// to find the overlapping kps in cd and FD i.e. same kps reconstructed twice in a event from FD and CD  
		if(select_efkpckp()){

			if(found_eFT && fkps.size() == 1 && ckps.size()==1){
				float fkpckp_deltap = fkps.get(0).p()-ckps.get(0).p();
				float fkpckp_deltatheta =  Math.toDegrees(fkps.get(0).theta())-Math.toDegrees(ckps.get(0).theta());
				float fkpckp_deltaphi = Math.toDegrees(fkps.get(0).phi())-Math.toDegrees(ckps.get(0).phi());
				hi_fkpckp_deltap.fill(fkpckp_deltap);
				hi_fkpckp_deltatheta.fill(fkpckp_deltatheta);
				hi_fkpckp_deltaphi.fill(fkpckp_deltaphi);
				if((fkpckp_deltatheta > -8.0) && (fkpckp_deltaphi > -25.0 || fkpckp_deltaphi < 10.0)){ //(fkpckp_deltap > -0.5 || fkpckp_deltap < 0.1) && 
					hi_fkpckp_deltap_p.fill(ckps.get(0).p(),fkpckp_deltap);///ckps.get(0).p()
					hi_fkpckp_deltap_theta.fill(Math.toDegrees(ckps.get(0).theta()),fkpckp_deltap);///ckps.get(0).p()
					hi_fkpckp_deltap_phi.fill(Math.toDegrees(ckps.get(0).phi()),fkpckp_deltap);///ckps.get(0).p()

					
					hi_fkpckp_deltatheta_p.fill(ckps.get(0).p(),fkpckp_deltatheta);///Math.toDegrees(ckps.get(0).theta())
					hi_fkpckp_deltatheta_theta.fill(Math.toDegrees(ckps.get(0).theta()),fkpckp_deltatheta);///Math.toDegrees(ckps.get(0).theta())
					hi_fkpckp_deltatheta_phi.fill(Math.toDegrees(ckps.get(0).phi()),fkpckp_deltatheta);///Math.toDegrees(ckps.get(0).theta())
						
					hi_fkpckp_deltaphi_p.fill(ckps.get(0).p(),fkpckp_deltaphi);///Math.toDegrees(ckps.get(0).phi())
					hi_fkpckp_deltaphi_theta.fill(Math.toDegrees(ckps.get(0).theta()),fkpckp_deltaphi);///Math.toDegrees(ckps.get(0).phi())
					hi_fkpckp_deltaphi_phi.fill(Math.toDegrees(ckps.get(0).phi()),fkpckp_deltaphi);///Math.toDegrees(ckps.get(0).phi())
					

					hi_fkpckp_deltap_withdpdtcut.fill(fkpckp_deltap);
					hi_fkpckp_deltatheta_withdpdtcut.fill(fkpckp_deltatheta);
					hi_fkpckp_deltaphi_withdpdtcut.fill(fkpckp_deltaphi);

				}
			}
		//	/*
			else if(found_eFT && kps.size() == 1 && fprots.size() == 1 && cprots.size()==1){
				float fkpckp_deltap = fprots.get(0).p()-cprots.get(0).p();
				float fkpckp_deltatheta =  Math.toDegrees(fprots.get(0).theta())-Math.toDegrees(cprots.get(0).theta());
				float fkpckp_deltaphi = Math.toDegrees(fprots.get(0).phi())-Math.toDegrees(cprots.get(0).phi());
				hi_fkpckp_deltap.fill(fkpckp_deltap);
				hi_fkpckp_deltatheta.fill(fkpckp_deltatheta);
				hi_fkpckp_deltaphi.fill(fkpckp_deltaphi);
				if((fkpckp_deltatheta > -8.0) && (fkpckp_deltaphi > -25.0 || fkpckp_deltaphi < 10.0)){ //(fkpckp_deltap > -0.5 || fkpckp_deltap < 0.1) && 
					hi_fkpckp_deltap_p.fill(cprots.get(0).p(),fkpckp_deltap);///cprots.get(0).p()
					hi_fkpckp_deltap_theta.fill(Math.toDegrees(cprots.get(0).theta()),fkpckp_deltap);///cprots.get(0).p()
					hi_fkpckp_deltap_phi.fill(Math.toDegrees(cprots.get(0).phi()),fkpckp_deltap);///cprots.get(0).p()

					
					hi_fkpckp_deltatheta_p.fill(cprots.get(0).p(),fkpckp_deltatheta);///Math.toDegrees(cprots.get(0).theta())
					hi_fkpckp_deltatheta_theta.fill(Math.toDegrees(cprots.get(0).theta()),fkpckp_deltatheta);///Math.toDegrees(cprots.get(0).theta())
					hi_fkpckp_deltatheta_phi.fill(Math.toDegrees(cprots.get(0).phi()),fkpckp_deltatheta);///Math.toDegrees(cprots.get(0).theta())
						
					hi_fkpckp_deltaphi_p.fill(cprots.get(0).p(),fkpckp_deltaphi);///Math.toDegrees(cprots.get(0).phi())
					hi_fkpckp_deltaphi_theta.fill(Math.toDegrees(cprots.get(0).theta()),fkpckp_deltaphi);///Math.toDegrees(cprots.get(0).phi())
					hi_fkpckp_deltaphi_phi.fill(Math.toDegrees(cprots.get(0).phi()),fkpckp_deltaphi);///Math.toDegrees(cprots.get(0).phi())
					
					hi_fkpckp_deltap_withdpdtcut.fill(fkpckp_deltap);
					hi_fkpckp_deltatheta_withdpdtcut.fill(fkpckp_deltatheta);
					hi_fkpckp_deltaphi_withdpdtcut.fill(fkpckp_deltaphi);

				}
			}

		 /*
			else if(found_eFT && kps.size() >= 1 && fpips.size() >= 1 && cpips.size()>=1){
				float fkpckp_deltap = fpips.get(0).p()-cpips.get(0).p();
				float fkpckp_deltatheta =  Math.toDegrees(fpips.get(0).theta())-Math.toDegrees(cpips.get(0).theta());
				float fkpckp_deltaphi = Math.toDegrees(fpips.get(0).phi())-Math.toDegrees(cpips.get(0).phi());
				hi_fkpckp_deltap.fill(fkpckp_deltap);
				hi_fkpckp_deltatheta.fill(fkpckp_deltatheta);
				hi_fkpckp_deltaphi.fill(fkpckp_deltaphi);

			}
			
			//*/
			
		}

		if(select_eckpckpckm()){
			hi_eckpckpckm_MM_eckpckp.fill(eckpckpckm_MM_eckpckp);
			hi_eckpckpckm_MM_eckpckpckm.fill(eckpckpckm_MM_eckpckpckm);
			hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm.fill(eckpckpckm_MM_eckpckpckm, eckpckpckm_MM_eckpckp);
			if(eckpckpckm_found_lambda){
					hi_eckpckpckm_MM_eckpckp_lam_evnt.fill(eckpckpckm_MM_eckpckp);
					hi_eckpckpckm_IM_kmlambda.fill(eckpckpckm_IM_kmlambda);
			}
		}
		
		if(select_ekpkp_req_km()) {
			//System.out.println(" kps.size() &&....... :) " + kps.size());
			//hi_ekpkpkm_MM_ekpkp.fill(ekpkpkm_MM_ekpkp);
			//hi_ekpkpkm_MM_ekpkpkm.fill(ekpkpkm_MM_ekpkpkm);

			//with electron common vertex correction
			hi_ekpkpkm_MM_ekpkp_nocorr.fill(ekpkpkm_MM_ekpkp_nocorr);
			hi_ekpkpkm_MM_ekpkpkm_nocorr.fill(ekpkpkm_MM_ekpkpkm_nocorr);

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

			hi_FT_e_t_f.fill(Math.toDegrees(Vftel_corrected.phi()), Math.toDegrees(Vftel_corrected.theta()));
			hi_FT_e_p_f.fill(Math.toDegrees(Vftel_corrected.phi()), Vftel_corrected.p());
			hi_FT_e_p_the.fill(Math.toDegrees(Vftel_corrected.theta()), Vftel_corrected.p());
			/*
			hi_rec_e_dp.fill((Vftel_corrected.p() - (float)mcels.get(0).p())/Vftel_corrected.p());
			hi_rec_e_dtheta.fill(Math.toDegrees(Vftel_corrected.theta()) - (float)Math.toDegrees(mcels.get(0).theta()));
			hi_rec_e_dphi.fill(Math.toDegrees(Vftel_corrected.phi()) - (float)Math.toDegrees(mcels.get(0).phi()));
			//dg_rec_electron.getH1F("hi_rec_e_dvx").fill(avg_vx - (float)mcels.get(0).vx());
			//dg_rec_electron.getH1F("hi_rec_e_dvy").fill(avg_vy - (float)mcels.get(0).vy());
			dg_rec_electron.getH1F("hi_rec_e_dvz").fill(avg_vz - (float)mcels.get(0).vz());
			dg_rec_e_resolution.getH2F("hi_rec_e_dp_p").fill(Vftel_corrected.p(), (Vftel_corrected.p() - (float)mcels.get(0).p())/Vftel_corrected.p());
			dg_rec_e_resolution.getH2F("hi_rec_e_dp_theta").fill(Math.toDegrees(Vftel_corrected.theta()), (Vftel_corrected.p() - (float)mcels.get(0).p())/Vftel_corrected.p());
			dg_rec_e_resolution.getH2F("hi_rec_e_dp_phi").fill(e_phi, (Vftel_corrected.p() - (float)mcels.get(0).p())/Vftel_corrected.p());
			dg_rec_e_resolution.getH2F("hi_rec_e_dp_vz").fill(avg_vz, (Vftel_corrected.p() - (float)mcels.get(0).p())/Vftel_corrected.p());
			dg_rec_e_resolution.getH2F("hi_rec_e_dtheta_p").fill(Vftel_corrected.p(), Math.toDegrees(Vftel_corrected.theta()) - (float)Math.toDegrees(mcels.get(0).theta()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dtheta_theta").fill(Math.toDegrees(Vftel_corrected.theta()), Math.toDegrees(Vftel_corrected.theta()) - (float)Math.toDegrees(mcels.get(0).theta()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dtheta_phi").fill(Math.toDegrees(Vftel_corrected.phi()), Math.toDegrees(Vftel_corrected.theta()) - (float)Math.toDegrees(mcels.get(0).theta()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dtheta_vz").fill(avg_vz, Math.toDegrees(Vftel_corrected.theta()) - (float)Math.toDegrees(mcels.get(0).theta()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dphi_p").fill(Vftel_corrected.p(), Math.toDegrees(Vftel_corrected.phi()) - (float)Math.toDegrees(mcels.get(0).phi()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dphi_theta").fill(Math.toDegrees(Vftel_corrected.theta()), Math.toDegrees(Vftel_corrected.phi()) - (float)Math.toDegrees(mcels.get(0).phi()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dphi_phi").fill(Math.toDegrees(Vftel_corrected.phi()), Math.toDegrees(Vftel_corrected.phi()) - (float)Math.toDegrees(mcels.get(0).phi()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dphi_vz").fill(avg_vz, Math.toDegrees(Vftel_corrected.phi()) - (float)Math.toDegrees(mcels.get(0).phi()));
			dg_rec_e_resolution.getH2F("hi_rec_e_dvz_p").fill(Vftel_corrected.p(), avg_vz - (float)mcels.get(0).vz());
			dg_rec_e_resolution.getH2F("hi_rec_e_dvz_theta").fill(Math.toDegrees(Vftel_corrected.theta()), avg_vz - (float)mcels.get(0).vz());
			dg_rec_e_resolution.getH2F("hi_rec_e_dvz_phi").fill(Math.toDegrees(Vftel_corrected.phi()), avg_vz - (float)mcels.get(0).vz());
			dg_rec_e_resolution.getH2F("hi_rec_e_dvz_vz").fill(avg_vz, avg_vz - (float)mcels.get(0).vz());
//*/

			
			//e resolution
			

			// kp1 resolution
			hi_rec_kp1_p_the.fill(rec_kp1_the, rec_kp1_p);
			hi_rec_kp1_p_phi.fill(rec_kp1_phi, rec_kp1_p);
			hi_rec_kp1_the_phi.fill(rec_kp1_phi, rec_kp1_the);
/*
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
//*/
			//kp2 resolution
			hi_rec_kp2_p_the.fill(rec_kp2_the, rec_kp2_p);
			hi_rec_kp2_p_phi.fill(rec_kp2_phi, rec_kp2_p);
			hi_rec_kp2_the_phi.fill(rec_kp2_phi, rec_kp2_the);
/*		
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
//*/
			//km resolution
			hi_rec_km_p_the.fill(rec_km_the, rec_km_p);
			hi_rec_km_p_phi.fill(rec_km_phi, rec_km_p);
			hi_rec_km_the_phi.fill(rec_km_phi, rec_km_the);
/*			
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
// */
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

		//myCanvas = new EmbeddedCanvasTabbed("REC-electron","electron-res","REC-kp1(fast)","kp1(fast-kp)-res","REC-kp2(slow)","kp2(slow-kp)-res", "REC-km","km-res","REC-Xi","MM-spectra", "FD-TOF", "CD-TOF", "VTime", "particle-Vz","TOF-t","TOF-path", "Counter", "Request");
		myCanvas = new EmbeddedCanvasTabbed("REC-electron","electron-res","REC-kp1(fast)","kp1(fast-kp)-res","REC-kp2(slow)","kp2(slow-kp)-res","REC-km","km-res","REC-Xi","MM-ekpkpkm","MM-ekpkpkm1","MM-scatter","MM-lam-evnt","IM-kmlambda","fkpckp","MM-spectra", "FD-TOF", "CD-TOF", "CD-Part", "VTime", "particle-Vz","TOF-t","TOF-path", "Counter", "Request");
		//myCanvas = new EmbeddedCanvas("REC-electron","electron-res","REC-kp1(fast)","kp1(fast-kp)-res","REC-kp2(slow)","kp2(slow-kp)-res", "REC-km","km-res","REC-Xi","MM-spectra", "FD-TOF", "CD-TOF", "VTime", "particle-Vz","TOF-t","TOF-path", "Counter", "Request");
		//reconstructed electron

		myCanvas.getCanvas("REC-electron").divide(3, 5);
		myCanvas.getCanvas("REC-electron").setSize(1600, 1000);
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

		/*
		EmbeddedCanvas rec_electron = new EmbeddedCanvas();
		rec_electron.setSize(1600, 1000);
		rec_electron.divide(3, 5);
		rec_electron.setGridX(false);
		rec_electron.setGridY(false);
		rec_electron.setAxisFontSize(18);
		rec_electron.setAxisTitleSize(24);
		rec_electron.draw(dg_rec_electron);
		rec_electron.getPad(0).getAxisZ().setLog(true);
		rec_electron.getPad(1).getAxisZ().setLog(true);
		rec_electron.getPad(2).getAxisZ().setLog(true);
		rec_electron.getPad(9).getAxisZ().setLog(true);
		rec_electron.getPad(12).getAxisZ().setLog(true);
		rec_electron.getPad(13).getAxisZ().setLog(true);
		rec_electron.getPad(14).getAxisZ().setLog(true);
		rec_electron.update();
		rec_electron.save("electron_canvas.png");
		*/

	//	/*
		myCanvas.getCanvas("electron-res").divide(4,4);
		myCanvas.getCanvas("electron-res").setSize(1600, 1000);
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
	//	*/
		//reconstructed kp (kp1;fast k+ and slow k+; kp2)

		myCanvas.getCanvas("REC-kp1(fast)").divide(3,4);
		myCanvas.getCanvas("REC-kp1(fast)").setSize(1600, 1000);
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

	//	/*
		// fast kaon resolution
		myCanvas.getCanvas("kp1(fast-kp)-res").divide(4,4);
		myCanvas.getCanvas("kp1(fast-kp)-res").setSize(1600, 1000);
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
	//	*/
		// slow kaon
		myCanvas.getCanvas("REC-kp2(slow)").divide(3,4);
		myCanvas.getCanvas("REC-kp2(slow)").setSize(1600, 1000);
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

	//	/*
		//slow kaon resolution
		myCanvas.getCanvas("kp2(slow-kp)-res").divide(4,4);
		myCanvas.getCanvas("kp2(slow-kp)-res").setSize(1600, 1000);
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
	
	//	*/
		//reconstructed km
		myCanvas.getCanvas("REC-km").divide(3,4);
		myCanvas.getCanvas("REC-km").setSize(1600, 1000);
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

	//	/*
		myCanvas.getCanvas("km-res").divide(4,4);
		myCanvas.getCanvas("km-res").setSize(1600, 1000);
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
	//	*/

		myCanvas.getCanvas("REC-Xi").divide(2,4);
		myCanvas.getCanvas("REC-Xi").setSize(1600, 1000);
		myCanvas.getCanvas("REC-Xi").setGridX(false);
		myCanvas.getCanvas("REC-Xi").setGridY(false);
		myCanvas.getCanvas("REC-Xi").setAxisFontSize(18);
		myCanvas.getCanvas("REC-Xi").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-Xi").draw(dg_rec_xi);
		//myCanvas.getCanvas("REC-Xi").getPad(1).draw(L_lamda);
		//myCanvas.getCanvas("REC-Xi").getPad(1).draw(L_sigma);

		myCanvas.getCanvas("MM-spectra").divide(3,3);
		myCanvas.getCanvas("MM-spectra").setSize(1600, 1000);
		myCanvas.getCanvas("MM-spectra").setGridX(false);
		myCanvas.getCanvas("MM-spectra").setGridY(false);
		myCanvas.getCanvas("MM-spectra").setAxisFontSize(18);
		myCanvas.getCanvas("MM-spectra").setAxisTitleSize(24);
		myCanvas.getCanvas("MM-spectra").draw(dg_mm);

		myCanvas.getCanvas("MM-ekpkpkm").divide(3, 2);
		myCanvas.getCanvas("MM-ekpkpkm").setSize(1600, 1000);
		myCanvas.getCanvas("MM-ekpkpkm").setGridX(false);
		myCanvas.getCanvas("MM-ekpkpkm").setGridY(false);
		myCanvas.getCanvas("MM-ekpkpkm").setAxisFontSize(18);
		myCanvas.getCanvas("MM-ekpkpkm").setAxisTitleSize(24);
		myCanvas.getCanvas("MM-ekpkpkm").draw(dg_mm_ekpkpkm);

		myCanvas.getCanvas("MM-ekpkpkm1").divide(3, 2);
		myCanvas.getCanvas("MM-ekpkpkm1").setSize(1600, 1000);
		myCanvas.getCanvas("MM-ekpkpkm1").setGridX(false);
		myCanvas.getCanvas("MM-ekpkpkm1").setGridY(false);
		myCanvas.getCanvas("MM-ekpkpkm1").setAxisFontSize(18);
		myCanvas.getCanvas("MM-ekpkpkm1").setAxisTitleSize(24);
		myCanvas.getCanvas("MM-ekpkpkm1").draw(dg_mm_ekpkpkm1);

		myCanvas.getCanvas("MM-scatter").divide(3, 2);
		myCanvas.getCanvas("MM-scatter").setSize(1600, 1000);
		myCanvas.getCanvas("MM-scatter").setGridX(false);
		myCanvas.getCanvas("MM-scatter").setGridY(false);
		myCanvas.getCanvas("MM-scatter").setAxisFontSize(18);
		myCanvas.getCanvas("MM-scatter").setAxisTitleSize(24);
		myCanvas.getCanvas("MM-scatter").draw(dg_mm_scatter);
		myCanvas.getCanvas("MM-scatter").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM-scatter").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM-scatter").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM-scatter").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM-scatter").getPad(4).getAxisZ().setLog(true);
		myCanvas.getCanvas("MM-scatter").getPad(5).getAxisZ().setLog(true);

	
		myCanvas.getCanvas("MM-lam-evnt").divide(3, 2);
		myCanvas.getCanvas("MM-lam-evnt").setSize(1600, 1000);
		myCanvas.getCanvas("MM-lam-evnt").setGridX(false);
		myCanvas.getCanvas("MM-lam-evnt").setGridY(false);
		myCanvas.getCanvas("MM-lam-evnt").setAxisFontSize(18);
		myCanvas.getCanvas("MM-lam-evnt").setAxisTitleSize(24);
		myCanvas.getCanvas("MM-lam-evnt").draw(dg_mm_ekpkpkm_lam_evnt);

		myCanvas.getCanvas("IM-kmlambda").divide(3, 2);
		myCanvas.getCanvas("IM-kmlambda").setSize(1600, 1000);
		myCanvas.getCanvas("IM-kmlambda").setGridX(false);
		myCanvas.getCanvas("IM-kmlambda").setGridY(false);
		myCanvas.getCanvas("IM-kmlambda").setAxisFontSize(18);
		myCanvas.getCanvas("IM-kmlambda").setAxisTitleSize(24);
		myCanvas.getCanvas("IM-kmlambda").draw(dg_m_kmlambda);

		myCanvas.getCanvas("fkpckp").divide(4,4);
		myCanvas.getCanvas("fkpckp").setSize(1600, 1000);
		myCanvas.getCanvas("fkpckp").setGridX(false);
		myCanvas.getCanvas("fkpckp").setGridY(false);
		myCanvas.getCanvas("fkpckp").setAxisFontSize(18);
		myCanvas.getCanvas("fkpckp").setAxisTitleSize(24);
		myCanvas.getCanvas("fkpckp").draw(dg_fkpckp_pthphi);
		//myCanvas.getCanvas("fkpckp").cd(5);
		//myCanvas.getCanvas("fkpckp").draw(meanX);

		myCanvas.getCanvas("FD-TOF").divide(5, 2);
		myCanvas.getCanvas("FD-TOF").setSize(1600, 1000);
		myCanvas.getCanvas("FD-TOF").setGridX(false);
		myCanvas.getCanvas("FD-TOF").setGridY(false);
		myCanvas.getCanvas("FD-TOF").setAxisFontSize(18);
		myCanvas.getCanvas("FD-TOF").setAxisTitleSize(24);
		myCanvas.getCanvas("FD-TOF").draw(dg_fdtof);
		myCanvas.getCanvas("FD-TOF").getPad(0).getAxisX().setRange(0.5, 5);
		myCanvas.getCanvas("FD-TOF").getPad(0).getAxisY().setRange(0.85, 1.05);
		myCanvas.getCanvas("FD-TOF").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(4).getAxisY().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(5).getAxisX().setRange(0.5, 5);
		myCanvas.getCanvas("FD-TOF").getPad(5).getAxisY().setRange(0.85, 1.05);
		myCanvas.getCanvas("FD-TOF").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("FD-TOF").getPad(9).getAxisY().setLog(true);
		//calculations from CD-TOF
		myCanvas.getCanvas("CD-TOF").divide(5, 2);
		myCanvas.getCanvas("CD-TOF").setSize(1600, 1000);
		myCanvas.getCanvas("CD-TOF").setGridX(false);
		myCanvas.getCanvas("CD-TOF").setGridY(false);
		myCanvas.getCanvas("CD-TOF").setAxisFontSize(18);
		myCanvas.getCanvas("CD-TOF").setAxisTitleSize(24);
		myCanvas.getCanvas("CD-TOF").draw(dg_cdtof);
		myCanvas.getCanvas("CD-TOF").getPad(0).getAxisX().setRange(0.2, 2.0);
		myCanvas.getCanvas("CD-TOF").getPad(0).getAxisY().setRange(0.1, 1.1);
		myCanvas.getCanvas("CD-TOF").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(2).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(3).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(4).getAxisY().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(5).getAxisX().setRange(0.2, 2.0);
		myCanvas.getCanvas("CD-TOF").getPad(5).getAxisY().setRange(0.1, 1.1);
		myCanvas.getCanvas("CD-TOF").getPad(5).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(6).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(7).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(8).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-TOF").getPad(9).getAxisY().setLog(true);

		myCanvas.getCanvas("CD-Part").divide(3, 3);
		myCanvas.getCanvas("CD-Part").setSize(1600, 1000);
		myCanvas.getCanvas("CD-Part").setGridX(false);
		myCanvas.getCanvas("CD-Part").setGridY(false);
		myCanvas.getCanvas("CD-Part").setAxisFontSize(18);
		myCanvas.getCanvas("CD-Part").setAxisTitleSize(24);
		myCanvas.getCanvas("CD-Part").draw(dg_cdPart);
		myCanvas.getCanvas("CD-Part").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-Part").getPad(4).getAxisZ().setLog(true);
		myCanvas.getCanvas("CD-Part").getPad(6).getAxisZ().setLog(true);

		// vertex time hadrons
		myCanvas.getCanvas("VTime").divide(5, 2);
		myCanvas.getCanvas("VTime").setSize(1600, 1000);
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
		myCanvas.getCanvas("particle-Vz").setSize(1600, 1000);
		myCanvas.getCanvas("particle-Vz").setGridX(false);
		myCanvas.getCanvas("particle-Vz").setGridY(false);
		myCanvas.getCanvas("particle-Vz").setAxisFontSize(18);
		myCanvas.getCanvas("particle-Vz").setAxisTitleSize(24);
		myCanvas.getCanvas("particle-Vz").draw(dg_vz);

		myCanvas.getCanvas("TOF-t").divide(5, 2);
		myCanvas.getCanvas("TOF-t").setSize(1600, 1000);
		myCanvas.getCanvas("TOF-t").setGridX(false);
		myCanvas.getCanvas("TOF-t").setGridY(false);
		myCanvas.getCanvas("TOF-t").setAxisFontSize(18);
		myCanvas.getCanvas("TOF-t").setAxisTitleSize(24);
		myCanvas.getCanvas("TOF-t").draw(dg_tof_t);

		myCanvas.getCanvas("TOF-path").divide(5, 2);
		myCanvas.getCanvas("TOF-path").setSize(1600, 1000);
		myCanvas.getCanvas("TOF-path").setGridX(false);
		myCanvas.getCanvas("TOF-path").setGridY(false);
		myCanvas.getCanvas("TOF-path").setAxisFontSize(18);
		myCanvas.getCanvas("TOF-path").setAxisTitleSize(24);
		myCanvas.getCanvas("TOF-path").draw(dg_tof_path);

		
		myCanvas.getCanvas("Counter").divide(5, 3);
		myCanvas.getCanvas("Counter").setSize(1600, 1000);
		myCanvas.getCanvas("Counter").setGridX(false);
		myCanvas.getCanvas("Counter").setGridY(false);
		myCanvas.getCanvas("Counter").setAxisFontSize(18);
		myCanvas.getCanvas("Counter").setAxisTitleSize(24);
		myCanvas.getCanvas("Counter").draw(dg_counter);


		myCanvas.getCanvas("Request").divide(2, 1);
		myCanvas.getCanvas("Request").setSize(1600, 1000);
		myCanvas.getCanvas("Request").setGridX(false);
		myCanvas.getCanvas("Request").setGridY(false);
		myCanvas.getCanvas("Request").setAxisFontSize(18);
		myCanvas.getCanvas("Request").setAxisTitleSize(24);
		myCanvas.getCanvas("Request").draw(dg_req);

		/*
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
		myCanvas.getCanvas("TOF-path").save("tof_pathi_overview.png");
		myCanvas.getCanvas("Counter").save("particle_multiplicity.png");
		myCanvas.getCanvas("Request").save("reconstructed_cascade_1820.png");
	*/

				
	}
	
	public void save() {
		
		//TDirectory dirout = new TDirectory();
		//dirout.mkdir("");
		//dirout.cd("/FTElec/");
		//dirout.addDataSet(hi_FT_e_t_f, hi_FT_e_p_f, hi_FT_W_Q2);
		myCanvas.getCanvas("REC-electron").save("ft_electron.png");
		//myCanvas.getCanvas("electron-res").save("ft_electron_resolution.png");
		myCanvas.getCanvas("REC-kp1(fast)").save("fast_kaon_plus.png");
		myCanvas.getCanvas("REC-kp2(slow)").save("slow_kaon_plus.png");
		//myCanvas.getCanvas("kp1(fast-kp)-res").save("fast_kaon_plus_resolution.png");
		//myCanvas.getCanvas("kp2(slow-kp)-res").save("slow_kaon_plus_resolution.png");
		myCanvas.getCanvas("REC-km").save("kaon_minus.png");
		//myCanvas.getCanvas("km-res").save("kaon_minus_resolution.png");
		myCanvas.getCanvas("REC-Xi").save("cascade.png");
		myCanvas.getCanvas("MM-ekpkpkm").save("ekpkpkm_mm_ekpkp_phasespace.png");
		myCanvas.getCanvas("MM-ekpkpkm1").save("ekpkpkm_mm_ekpkpkm_phasespace.png");
		myCanvas.getCanvas("MM-scatter").save("ekpkpkm_mm_ekpkp_mm_ekpkpkm_scatter.png");
		myCanvas.getCanvas("MM-lam-evnt").save("ekpkpkm_mm_ekpkp_lambdaevent.png");
		myCanvas.getCanvas("IM-kmlambda").save("ekpkpkm_im_kmlambda.png")
		myCanvas.getCanvas("fkpckp").save("detector_overlap_fkpckp.png");
		myCanvas.getCanvas("MM-spectra").save("mm_ekpkpkm.png");
		myCanvas.getCanvas("CD-TOF").save("cd_beta_p_overview.png");
		myCanvas.getCanvas("FD-TOF").save("fd_beta_p_overview.png");
		myCanvas.getCanvas("VTime").save("particle_dtime_overview.png");
		myCanvas.getCanvas("particle-Vz").save("particle_vz_overview.png");
		myCanvas.getCanvas("TOF-t").save("tof_time_overview.png");
		myCanvas.getCanvas("TOF-path").save("tof_pathi_overview.png");
		myCanvas.getCanvas("Counter").save("particle_multiplicity.png");
		myCanvas.getCanvas("Request").save("reconstructed_cascade_1820.png");
		

	}

	public void write() {
		
		TDirectory dirout = new TDirectory();
		dirout.mkdir("/plots");
		dirout.cd("/plots");
		dirout.addDataSet(hi_FT_e_t_f,hi_FT_e_p_f,hi_FT_e_p_the,hi_FT_Q2,hi_FT_W,hi_FT_W_Q2,hi_virphoton,hi_rec_kp1_p_the,hi_rec_kp1_p_phi,hi_rec_kp1_the_phi,hi_rec_kp2_p_the,hi_rec_kp2_p_phi,hi_rec_kp2_the_phi,hi_rec_km_p_the,hi_rec_km_p_phi,hi_rec_km_the_phi);
		//dirout.save("test.hipo");
		dirout.writeFile("histFiles.hipo");

		//TDirectory rdir = new TDirectory();
        //rdir.load("out_FT.hipo");
        //rdir.ls();
        //TBrowser t = new TBrowser(rdir);

	}

	public void showplots() {

		JFrame frame = new JFrame("RGA-Fall2018");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(1600, 1000);
		//frame.setSize(1800, 1200);
		frame.add(myCanvas);
		frame.setLocationRelativeTo(null);
		//frame.pack();
		//frame.setMinimumSize(new Dimension(300,300))
		frame.setVisible(true);
		//frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		//frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		//frame.setDefaultCloseOperation();

	}

	static void eventMixing(ArrayList<LorentzVector> els, ArrayList<LorentzVector> fast_kps, ArrayList<LorentzVector> slow_kps, ArrayList<LorentzVector> kms){

		//System.out.println("Mp "+ Mp +" el-size: "+els.size()+ " fast-kp-size: "+fast_kps.size()+ " slow-kp-size: "+ slow_kps.size()+ " km-size: "+ kms.size());
		
		///*
		//int i = 0;
		//while(i<kms.size()){
		for(int i = 0; i < kms.size(); i++){
			LorentzVector lv_gs = new LorentzVector();
			lv_gs.add(VB);lv_gs.add(VT);
			for(int j = 0; j < fast_kps.size(); j++){
				if (i == j) {continue;}
				LorentzVector lv_mm_fkpskp = new LorentzVector();lv_mm_fkpskp.add(lv_gs);lv_mm_fkpskp.sub(els.get(i));lv_mm_fkpskp.sub(fast_kps.get(j));lv_mm_fkpskp.sub(slow_kps.get(j));
				double mm_fkpskp_mixed = lv_mm_fkpskp.mass();
				LorentzVector lv_mm_fkpskpkm_mixed = new LorentzVector();lv_mm_fkpskpkm_mixed.add(lv_mm_fkpskp);lv_mm_fkpskpkm_mixed.sub(kms.get(i));
				double mm_fkpskpkm_mixed = lv_mm_fkpskpkm_mixed.mass();
				me_outFile.append(mm_fkpskpkm_mixed+" "+mm_fkpskp_mixed+"\n");
				//System.out.println("MM(ekpkpkm)"+ mm_fkpskpkm_mixed +" Mp "+ Mp +" el-size: "+els.size()+ " fast-kp-size: "+fast_kps.size()+ " slow-kp-size: "+ slow_kps.size()+ " km-size: "+ kms.size());
			}
		//i++; while loop
		}
		
		//*/
	}

	static void threeKaon_eventMixing(ArrayList<LorentzVector> els, ArrayList<LorentzVector> fast_kps, ArrayList<LorentzVector> slow_kps, ArrayList<LorentzVector> kms){

		System.out.println("Mp "+ Mp +" el-size: "+els.size()+ " fast-kp-size: "+fast_kps.size()+ " slow-kp-size: "+ slow_kps.size()+ " km-size: "+ kms.size());
		
		///*
		for(int i = 0; i < kms.size(); i++){
			LorentzVector lv_gs = new LorentzVector();
			lv_gs.add(VB);lv_gs.add(VT);
			for(int j = 0; j < fast_kps.size(); j++){
				for(int k = 0; k < slow_kps.size(); k++){
					if (i == j && j == k) {continue;}
					LorentzVector lv_mm_fkpskp = new LorentzVector();lv_mm_fkpskp.add(lv_gs);lv_mm_fkpskp.sub(els.get(i));lv_mm_fkpskp.sub(fast_kps.get(j));lv_mm_fkpskp.sub(slow_kps.get(k));
					double mm_fkpskp_mixed = lv_mm_fkpskp.mass();
					LorentzVector lv_mm_fkpskpkm_mixed = new LorentzVector();lv_mm_fkpskpkm_mixed.add(lv_mm_fkpskp);lv_mm_fkpskpkm_mixed.sub(kms.get(i));
					double mm_fkpskpkm_mixed = lv_mm_fkpskpkm_mixed.mass();
					me_outFile.append(mm_fkpskpkm_mixed+" "+mm_fkpskp_mixed+"\n");
					//System.out.println("MM(ekpkpkm)"+ mm_fkpskpkm_mixed +" Mp "+ Mp +" el-size: "+els.size()+ " fast-kp-size: "+fast_kps.size()+ " slow-kp-size: "+ slow_kps.size()+ " km-size: "+ kms.size());

				}
				
			}
		}
		
		//*/
	}
	
	static void twoKaon_eventMixing(ArrayList<LorentzVector> els, ArrayList<LorentzVector> fast_kps, ArrayList<LorentzVector> slow_kps){
		for(int i = 0; i < fast_kps.size(); i++){
			LorentzVector lv_gs = new LorentzVector();
			lv_gs.add(VB);lv_gs.add(VT);
			for(int j = 0; j < slow_kps.size(); j++){
				if (i == j) {continue;}
				LorentzVector lv_mm_fkpskp_mixed = new LorentzVector();
				lv_mm_fkpskp_mixed.add(lv_gs);lv_mm_fkpskp_mixed.sub(els.get(i));lv_mm_fkpskp_mixed.sub(fast_kps.get(i));lv_mm_fkpskp_mixed.sub(slow_kps.get(j));
				double mm_fkpskp_mixed = lv_mm_fkpskp_mixed.mass();
				me_twoKaon_outFile.append(mm_fkpskp_mixed+"\n");
			}
		}
	}
		
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	public static void main(String[] args) {
		//System.setProperty("java.awt.headless", "true"); 
		GStyle.setPalette("kRainBow");
		GStyle.getH1FAttributes().setOptStat("1110");
        //GStyle.getFunctionAttributes().setOptStat("1100");

		int count = 0;
		//long maxevents = 8000;
		int maxevents = 1000000/2;
		//long maxevents = 200000000000;
		pionsaskaons_ekpkpreqkmavgvz_data_v1 ana = new pionsaskaons_ekpkpreqkmavgvz_data_v1();
		System.out.println(String.format(">>> files from list %s >>>", args[0]));
		String filelist;// = "list_of_files.txt";
		filelist = args[0];
		runType = args[1];
		ArrayList<String> toProcessFileNames = new ArrayList<String>();//List<String> toProcessFileNames = new ArrayList<String>();
		File file = new File(filelist);
		Scanner read;
        try {
                read = new Scanner(file);
                do {

                		if (runType == "data"){

                			String filename = read.next();
                            toProcessFileNames.add(filename);

                		} else if (runType == "mc"){

                			String filename = read.next();
                			//toProcessFileNames.add(filename + "/dst.hipo"); //for reading mc files
                			toProcessFileNames.add(filename);

                		} else {
                			println("You provided wrong argument for runType. Please specify mc or data")
                		}
                	

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

        //HipoDataSync  fileWriter1 = new HipoDataSync();
        //fileWriter1.open("filter_efkpfkpfkm.hipo");

         /*
        HipoDataSync  fileWriter1 = new HipoDataSync();
        HipoDataSync  fileWriter2 = new HipoDataSync();
        HipoDataSync  fileWriter3 = new HipoDataSync();
        HipoDataSync  fileWriter4 = new HipoDataSync();
        HipoDataSync  fileWriter5 = new HipoDataSync();
        HipoDataSync  fileWriter6 = new HipoDataSync();
        fileWriter1.open("filter_efkpfkpfkm.hipo");
        fileWriter2.open("filter_efkpckpfkm.hipo");
        fileWriter3.open("filter_eckpckpfkm.hipo");
        fileWriter4.open("filter_efkpfkpckm.hipo");
        fileWriter5.open("filter_efkpckpckm.hipo");
        fileWriter6.open("filter_eckpckpckm.hipo");
    //    */

       /*
       	HipoDataSync  fileWriter1 = new HipoDataSync();
       	HipoDataSync  fileWriter2 = new HipoDataSync();
       	HipoDataSync  fileWriter3 = new HipoDataSync();
		fileWriter1.open("filter_fkpfkpfkm.hipo");
		fileWriter2.open("filter_fkpckpfkm.hipo");
		fileWriter3.open("filter_fkpckp.hipo");
		//  */


		mcoutFile = new File("gen_output_file.txt");
       	mcoutFile.bytes = new byte[0]

       	outFile = new File("output_file.txt");
       	outFile.bytes = new byte[0]

       	me_outFile = new File("me_output_file.txt");
       	me_outFile.bytes = new byte[0]

       	me_twoKaon_outFile = new File("me_twoKaon_file.txt");
       	me_twoKaon_outFile.bytes = new byte[0]

       	//public static ArrayList<LorentzVector> lv_my_ftels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
       	lv_my_ftels = new ArrayList<LorentzVector>();
		lv_my_fast_fkps = new ArrayList<LorentzVector>();
		lv_my_slow_fkps = new ArrayList<LorentzVector>();
		lv_my_fkms = new ArrayList<LorentzVector>();

		//public static ArrayList<LorentzVector> lv_my_xi_ftels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms;
		lv_my_xi_ftels = new ArrayList<LorentzVector>();
		lv_my_xi_fast_fkps = new ArrayList<LorentzVector>();
		lv_my_xi_slow_fkps = new ArrayList<LorentzVector>();
		lv_my_xi_fkms = new ArrayList<LorentzVector>();

        for (String runstrg : toProcessFileNames) if(count < maxevents){ //if(count < maxevents)
        	progresscount++;
            System.out.println(String.format(">>>> file %s >>>>", runstrg));
            File varTmpDir = new File(runstrg);
            if(!varTmpDir.exists()){System.out.println("FILE DOES NOT EXIST");continue;}
            System.out.println("READING NOW " + runstrg);

            if(runstrg == "/Users/akhanal/work/Expt_phys/simulation/xi_1820/rga_spring2019_in_skim11/plot/filter_efkpfkpfkm.hipo" || runstrg == "/Users/akhanal/work/Expt_phys/simulation/xi_1820/rga_spring2019_in_skim11/plot/filter_efkpckpfkm.hipo"){
            	System.out.println("SPRING 2019 FILE DOES EXIST SETTING Eb = 10.1998 GeV");
            	Eb = 10.1998f;
            	VB = new LorentzVector(0, 0, Eb, Eb);
            } else {
            	//Eb = 10.1998f;
            	//Eb = 10.6f;
            	Eb = 10.604f;
            	VB = new LorentzVector(0, 0, Eb, Eb);
            }

            HipoDataSource reader = new HipoDataSource();
            reader.open(runstrg);
            int filecount = 0;
            while(reader.hasEvent() && count < maxevents) { //while(reader.hasEvent()) && count < maxevents

            	//DataEvent event = reader.getNextEvent();
            	HipoDataEvent event = reader.getNextEvent();


            	runnum = event.getBank("RUN::config").getInt('run',0);
            	evnum = event.getBank("RUN::config").getInt('event',0);

            	if (runType == "mc"){
					//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
					Random rand = new Random();
					// fix smearing factor for FD Kps to 0.35% derived using events with electron in FD
    				smearFactor = 0.3*rand.nextGaussian()/100; 
				}

            	//Timothy's analysis fitter
            	//GenericKinematicFitter research_fitter = new analysis_fitter(Eb);
            	//research_fitter = new analysis_fitter(Eb); //research fitter defined as public static GenericKinematicFitter to use Timothys functions
            	//PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);

            	//if(research_fitter.banks_test(event)){
            	//	ana.processEvent(event);
            	//}


                ana.processEvent(event);


                //if(ana.select_efkpfkpfkm()){
                //	fileWriter1.writeEvent(event);
                //}
                /*
                if(ana.select_efkpfkpfkm()){
                	fileWriter1.writeEvent(event);
                }
                if(ana.select_efkpckpfkm()){
                	fileWriter2.writeEvent(event);
                }
                if(ana.select_eckpckpfkm()){
                	fileWriter3.writeEvent(event);
                }
                if(ana.select_efkpfkpckm()){
                	fileWriter4.writeEvent(event);
                }
                if(ana.select_efkpckpckm()){
                	fileWriter5.writeEvent(event);
                }
                if(ana.select_eckpckpckm()){
                	fileWriter6.writeEvent(event);
                }

			//	*/

               /*
                
                if(ana.select_efkpckp()){
                	fileWriter3.writeEvent(event);
                }
				//*/

                filecount++;count++;
                //for()

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

        //public static ArrayList<LorentzVector> lv_my_ftels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
        //eventMixing(lv_my_ftels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms);

        //public static ArrayList<LorentzVector> lv_my_xi_ftels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms;
        //twoKaon_eventMixing(lv_my_xi_ftels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps);


        //threeKaon_eventMixing(lv_my_fdels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms);
        threeKaon_eventMixing(lv_my_xi_ftels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms);

        //fileWriter1.close();
        /*
        fileWriter1.close();
		fileWriter2.close();
		fileWriter3.close();
		fileWriter4.close();
		fileWriter5.close();
		fileWriter6.close();
	//	*/
        /*
        fileWriter1.close();
        fileWriter2.close();
        fileWriter3.close();
		//*/

		System.out.println("Total events : " + count);
		//ana.analyze();
		/*
		ana.plot();
		ana.plotGraph();
		ana.showplots();
	//	ana.save();
		//*/
		//ana.save();
		
		// for ifarm
	//	/*
		ana.plot();
		ana.plotGraph();
		ana.showplots();
		//ana.save();
		//*/
		
		//ana.write();
		System.out.println("Good Bye !!!!!");
	}
		
	
	
}
