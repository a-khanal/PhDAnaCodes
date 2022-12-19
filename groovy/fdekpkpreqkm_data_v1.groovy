import java.io.*;
import java.util.*;
import org.jlab.io.hipo.*;
import java.text.SimpleDateFormat;
import org.jlab.jnp.physics.*;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.pdg.PDGParticle;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.physics.PhysicsEvent;
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
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import javax.swing.JFrame;
//import org.jlab.jroot.ROOTFile;
import org.jlab.clas.physics.GenericKinematicFitter;
import extended_kinematic_fitters.*;

/**
 * @author akhanal
 *
 */


public class fdekpkpreqkm_data_v1 {
	//public int NFTElec;
	


	public static float  Eb, Mp;
	public static LorentzVector VB;
	public static LorentzVector VT;

	//Timothy's analysis fitter
	public static GenericKinematicFitter research_fitter;
	public static File outFile, me_outFile, me_twoKaon_outFile;
	public static File mcoutFile;
	public static String runType;

	public static ArrayList<LorentzVector> lv_my_fdels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
	public static ArrayList<LorentzVector> lv_my_xi_fdels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms;

	//public float EB, Eb, Mp;
	public float STT, RFT, FTSTT, vt;
	public static int runnum, evnum;
	public static double smearFactor;
	
	public Particle Vprot, Vpip, Vpim, Vkp, Vkm;
	public Particle Vprotc, Vpipc, Vpimc, Vkpc, Vkmc;
	public Particle Vftel_corrected;
	public LorentzVector Ve_consCorn;
	public LorentzVector Ve, VGS, VhadronSystm, Vpim_correct, Vkp_correct;
	public Particle Vfde;

	public Vector3D e_ftCal_hitPosition;
	
	public boolean found_fastpip, found_slowpip;
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
	public F1D fn_fkp_deltat;
	
	
	public DataLine L_lamda, L_sigma;
	
	
	public H1F H_FT_W, H_FT_Q2, H_virphoton;
	public H1F  hi_ekpkp_MM_req_pim, hi_efkpckp_MM_req_pim, hi_efkpckp_MM_req_cpim, hi_efkpckp_MM_req_km, hi_efkpckp_MM_req_ckm;
	// particles vz, path_FTOF1b, time_FTOF1b
	
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

	public H1F hi_fd_pos_mass, hi_fd_neg_mass, hi_cd_pos_mass, hi_cd_neg_mass;
	public H2F hi_FD_pos_beta_mom, hi_FD_neg_beta_mom, hi_FD_neutral_beta_mom;
	public H2F hi_FD_pos_mass_mom, hi_FD_neg_mass_mom, hi_FD_neutral_mass_mom;
	public H2F hi_FD_pos_mass_the, hi_FD_neg_mass_the, hi_FD_neutral_mass_the;
	public H2F hi_FD_pos_mass_phi, hi_FD_neg_mass_phi, hi_FD_neutral_mass_phi;
	public H2F hi_CD_pos_beta_mom, hi_CD_neg_beta_mom, hi_CD_neutral_beta_mom;
	public H2F hi_CD_pos_mass_mom, hi_CD_neg_mass_mom, hi_CD_neutral_mass_mom;
	public H2F hi_CD_pos_mass_the, hi_CD_neg_mass_the, hi_CD_neutral_mass_the;
	public H2F hi_CD_pos_mass_phi, hi_CD_neg_mass_phi, hi_CD_neutral_mass_phi;

	public ArrayList<Particle> fpips, fpims, fkps, fkms, fprots; // = new ArrayList<Particle>();
	public ArrayList<Particle> cpips, cpims, ckps, ckms, cprots;
	public ArrayList<Particle> fdels;
	public ArrayList<Particle> ftels, pips, pims, kps, kms, prots;
	public ArrayList<Particle> mcels, mckps, mckms;

	// for e pip pip pim detected

	public float epippippim_mm_epippippim, epippippim_m_slowpippim, epippippim_m_fastpippim;
	public H1F hi_epippippim_mm_epippippim, hi_epippippim_m_slowpippim, hi_epippippim_m_fastpippim;


    public Particle mckp1;
    public Particle mckp2;
	// for e kp kp km detected
	public Particle reckp1;
    public Particle reckp2;
	public float rec_kp_p, rec_kp_the, rec_kp_phi, rec_kp_vz, rec_prot_p, rec_prot_the, rec_prot_phi, rec_prot_vz;
	public float rec_km_p, rec_km_the, rec_km_phi, rec_km_vz;

	public float efkpfkpfkm_MM_efkpfkp,  efkpckpfkm_MM_efkpckp,  eckpckpfkm_MM_eckpckp;  
	public float efkpfkpckm_MM_efkpfkp,  efkpckpckm_MM_efkpckp,  eckpckpckm_MM_eckpckp;
	public float efkpfkpfkm_MM_efkpfkpfkm,  efkpckpfkm_MM_efkpckpfkm, efkpckpfkm_MM_efkpckpfkm_nocorn, eckpckpfkm_MM_eckpckpfkm;  
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

    public F1D fn_lamsig_fit, fn_lam_fit, fn_sig_fit, fn_bg_fit; // in the MM(eK+K+K-) spectrum

    public F1D fn2_deltaphi_phi1, fn2_deltaphi_phi2, fn2_deltaphi_phi3;
    public F1D fn2_deltatheta_theta, fn2_deltatheta_p;
	
	public F1D fn_rec_e_dp, fn_rec_e_dtheta, fn_rec_e_dphi, fn_rec_e_dvx, fn_rec_e_dvy, fn_rec_e_dvz;
	public F1D fn_rec_kp_dp, fn_rec_kp_dtheta, fn_rec_kp_dphi, fn_rec_kp_dvx, fn_rec_kp_dvy, fn_rec_kp_dvz;
	public F1D fn_rec_prot_dp, fn_rec_prot_dtheta, fn_rec_prot_dphi, fn_rec_prot_dvx, fn_rec_prot_dvy, fn_rec_prot_dvz;
	public F1D fn_rec_km_dp, fn_rec_km_dtheta, fn_rec_km_dphi, fn_rec_km_dvx, fn_rec_km_dvy, fn_rec_km_dvz;

	public H2F hi_FT_e_beta_mom;
	public H2F H_FT_e_t_f, H_FT_e_p_f, H_FT_e_p_the;
	public H1F hi_rec_e_dp, hi_rec_e_dtheta, hi_rec_e_dphi;
	public H2F hi_rec_e_dp_p, hi_rec_e_dp_theta, hi_rec_e_dp_phi, hi_rec_e_dp_vz, hi_rec_e_dtheta_p, hi_rec_e_dtheta_theta, hi_rec_e_dtheta_phi, hi_rec_e_dtheta_vz, hi_rec_e_dphi_p, hi_rec_e_dphi_theta, hi_rec_e_dphi_phi, hi_rec_e_dphi_vz, hi_rec_e_dvz_p, hi_rec_e_dvz_theta, hi_rec_e_dvz_phi, hi_rec_e_dvz_vz;
	public H2F H_FT_W_Q2, H_FT_e_xB_Q2;
	public H2F hi_W_cd_pro_the;

	// momentum, theta and phi resolution for kp1, kp2 and km
	public H1F hi_rec_kp_dp, hi_rec_kp_dtheta, hi_rec_kp_dphi;
	public H2F hi_rec_kp_dp_p, hi_rec_kp_dp_theta, hi_rec_kp_dp_phi, hi_rec_kp_dp_vz, hi_rec_kp_dtheta_p, hi_rec_kp_dtheta_theta, hi_rec_kp_dtheta_phi, hi_rec_kp_dtheta_vz, hi_rec_kp_dphi_p, hi_rec_kp_dphi_theta, hi_rec_kp_dphi_phi, hi_rec_kp_dphi_vz, hi_rec_kp_dvz_p, hi_rec_kp_dvz_theta, hi_rec_kp_dvz_phi, hi_rec_kp_dvz_vz;
	public H1F hi_rec_prot_dp, hi_rec_prot_dtheta, hi_rec_prot_dphi;
	public H2F hi_rec_prot_dp_p, hi_rec_prot_dp_theta, hi_rec_prot_dp_phi, hi_rec_prot_dp_vz, hi_rec_prot_dtheta_p, hi_rec_prot_dtheta_theta, hi_rec_prot_dtheta_phi, hi_rec_prot_dtheta_vz, hi_rec_prot_dphi_p, hi_rec_prot_dphi_theta, hi_rec_prot_dphi_phi, hi_rec_prot_dphi_vz, hi_rec_prot_dvz_p, hi_rec_prot_dvz_theta, hi_rec_prot_dvz_phi, hi_rec_prot_dvz_vz;
	public H1F hi_rec_km_dp, hi_rec_km_dtheta, hi_rec_km_dphi;
	public H2F hi_rec_km_dp_p, hi_rec_km_dp_theta, hi_rec_km_dp_phi, hi_rec_km_dp_vz, hi_rec_km_dtheta_p, hi_rec_km_dtheta_theta, hi_rec_km_dtheta_phi, hi_rec_km_dtheta_vz, hi_rec_km_dphi_p, hi_rec_km_dphi_theta, hi_rec_km_dphi_phi, hi_rec_km_dphi_vz, hi_rec_km_dvz_p, hi_rec_km_dvz_theta, hi_rec_km_dvz_phi, hi_rec_km_dvz_vz;


	// scatter plotes for hadrons 
	public H2F hi_rec_kp1_p_the, hi_rec_kp2_p_the, hi_rec_km_p_the, hi_mc_e_p_the, hi_mc_kp1_p_the, hi_mc_kp2_p_the, hi_mc_km_p_the;
	public H2F hi_rec_kp1_p_phi, hi_rec_kp2_p_phi, hi_rec_km_p_phi, hi_mc_e_p_phi, hi_mc_kp1_p_phi, hi_mc_kp2_p_phi, hi_mc_km_p_phi;
	public H2F hi_rec_kp1_the_phi, hi_rec_kp2_the_phi, hi_rec_km_the_phi, hi_mc_e_the_phi, hi_mc_kp1_the_phi, hi_mc_kp2_the_phi, hi_mc_km_the_phi;

	public H1F hi_rec_kps_dtheta;

	public H2F hi_rec_kps_dtheta_slowkp_theta, hi_rec_kps_dtheta_fastkp_theta, hi_rec_kps_dtheta_slowkp_phi, hi_rec_kps_dtheta_fastkp_phi, hi_rec_kps_dtheta_slowkp_p, hi_rec_kps_dtheta_fastkp_p;


	public H1F hi_rec_fde_vz, hi_rec_kp1_vz, hi_rec_kp2_vz, hi_rec_km_vz;

	public H2F hi_rec_fde_theta_vz, hi_rec_kp1_theta_vz, hi_rec_kp2_theta_vz, hi_rec_km_theta_vz;


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

	public H2F hi_efkpfkpfkm_MM_efastkpaspipslowkpfkm, hi_efkpfkpfkm_MM_eslowkpaspipfastkpfkm;// considering fast kp as pip;

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

	public EmbeddedCanvas rec_electron;
	public EmbeddedCanvasTabbed myCanvas;
	//public EmbeddedCanvas myCanvas;
	public DataGroup dg_kinematics, dg_rec_electron, dg_rec_kp1, dg_rec_kp2, dg_rec_kps, dg_rec_km, dg_rec_xi, dg_rec_p, dg_rec_pim, dg_vtime, dg_fdtof, dg_cdtof, dg_cdPart, dg_vz, dg_tof_t, dg_tof_path, dg_counter;
	public DataGroup dg_rec_e_resolution, dg_rec_kp1_resolution, dg_rec_kp2_resolution, dg_rec_km_resolution;
	public DataGroup dg_mm, dg_mm_ekpkpkm, dg_mm_ekpkpkm1, dg_rec_kpaspip, dg_mm_scatter, dg_mm_ekpkpkm_lam_evnt, dg_m_kmlambda, dg_fkpckp_pthphi;
	public DataGroup dg_req;


	
	
	public fdekpkpreqkm_data_v1() {

		final int RED = 2;
		final int BLUE = 9;
		final int LIGHTGREEN = 3;
		final int LIGHTBROWN = 45;
		final int PINK = 46;
	//	NFTElec = 0;
	//	Eb = 10.575f;
	//	Eb = 10.604f; //RGA fall2018 beam energy
	//	Eb = 7.54626f;
	//	Eb = 6.535f
	//	Eb = 10.1998; //RGA spring2019 beam energy
		Mp = (float) PDGDatabase.getParticleById(2212).mass();
	//	Mp = 0.93827f;
		
	//	VB = new LorentzVector(0, 0, Eb, Eb);
		VT = new LorentzVector(0, 0, 0, Mp);

		
		// theoretical 1D functions for proton, kaon and pion
		
		//F_prot_beta_mom = new F1D("F-prot-beta-mom", "x/sqrt(0.93827*0.93827+x*x)", 0.3, 4.0);
		F_prot_beta_mom = new F1D("F-prot-beta-mom", "x/sqrt(0.93827*0.93827+x*x)", 0.1, 5.0);
		//F_prot_beta_mom.setRange(0.85, 1.05);
		F_prot_beta_mom.setLineWidth(2);
		F_prot_beta_mom.setLineColor(BLUE);
		//F_kp_beta_mom = new F1D("F-kp-beta-mom", "x/sqrt(0.49367*0.49367+x*x)", 0.3, 4.0);
		F_kp_beta_mom = new F1D("F-kp-beta-mom", "x/sqrt(0.49367*0.49367+x*x)", 0.1, 5.0);
		F_kp_beta_mom.setLineColor(BLUE);
		F_kp_beta_mom.setLineWidth(2);
		//F_pip_beta_mom = new F1D("F-pip-beta-mom", "x/sqrt(0.13957*0.13957+x*x)", 0.3, 4.0);
		F_pip_beta_mom = new F1D("F-pip-beta-mom", "x/sqrt(0.13957*0.13957+x*x)", 0.1, 5.0);
		F_pip_beta_mom.setLineColor(BLUE);
		F_pip_beta_mom.setLineWidth(2);

		/*
		// fit functions for Lambda, Sigma0, and Bg from event mixing in MM(eK+K+K-)
		//fn_lam_fit = new F1D("fn-lam-fit", "[amp]*gaus(x,[mean],[sigma])", 0, 5);
		fn_lam_fit = new F1D("fn-lam-fit", "[a]*gaus(x,[m],[s])", 0, 5);
		//fn_lam_fit.setParameters(29.97, 1.118, 0.03481) // for f18outfte data set
		fn_lam_fit.setParameters(31.15, 1.129, 0.02016) // for f18outfde data set

		//fn_sig_fit = new F1D("fn-sig-fit", "[amp]*gaus(x,[mean],[sigma])", 0, 5);
		fn_sig_fit = new F1D("fn-sig-fit", "[a]*gaus(x,[m],[s])", 0, 5);
		//fn_sig_fit.setParameters(13.81, 1.197, 0.03000) // for f18outfte data set
		fn_sig_fit.setParameters(15.00, 1.200, 0.01986) // for f18outfde data set

		fn_bg_fit = new F1D("fn-bg-fit", "([a]+[b]*x+[c]*x*x+[d]*x*x*x)*[e]", 0, 5);
		//fn_bg_fit.setParameters(65463, -211643, 260128, -76747, 0.0003136); // for f18outfte data set
		fn_bg_fit.setParameters(2323, -7075, 16887, -6546, 0.0009359); // for f18outfde data set
		//*/

		// fit double gaussian + bg (poly)
		fn_lamsig_fit = new F1D("fn-lamsig-fit", "[a1]*gaus(x,[m1],[s1]) + [a2]*gaus(x,[m2],[s2])", 0, 5);
		//fn_lamsig_fit.setParameters(31.1642, 1.12943, 0.0207800, 15.0000, 1.19955, 0.0182951);// for f18outfde data set
		//fn_lamsig_fit.setParameters(23.0000, 1.12782, 0.0255499, 4.00000, 1.20123, 0.0303990);// for f18infde data set
		fn_lamsig_fit.setParameters(24.0000, 1.12321, 0.0205438, 11.7659, 1.18969, 0.0275730);// for s19infde data set
		fn_lam_fit = new F1D("fn-lam-fit", "[a]*gaus(x,[m],[s])", 0, 5);
		//fn_lam_fit.setParameters(31.1642, 1.12943, 0.0207800); // for f18outfde data set
		//fn_lam_fit.setParameters(23.0000, 1.12782, 0.0255499); // for f18infde data set
		fn_lam_fit.setParameters(24.0000, 1.12321, 0.0205438); // for s19infde data set
		fn_sig_fit = new F1D("fn-sig-fit", "[a]*gaus(x,[m],[s])", 0, 5);
		//fn_sig_fit.setParameters(15.0000, 1.19955, 0.0182951) // for f18outfde data set
		//fn_sig_fit.setParameters(4.00000, 1.20123, 0.0303990); // for f18infde data set
		fn_sig_fit.setParameters(11.7659, 1.18969, 0.0275730); // for s19infde data set
		fn_bg_fit = new F1D("fn-bg-fit", "[a]+[b]*x", 0, 5);
		//fn_bg_fit.setParameters(-8.48096, 12.6331); // for f18outfde data set
		//fn_bg_fit.setParameters(-0.995676, 3.56739); // for f18infte data set
		fn_bg_fit.setParameters(-6.60240, 8.82551); // for s19infte data set



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
		

	// FD electron overview

		//reconstructed
		H_FT_e_t_f = new H2F("H-FT-e-t-f", "H-FT-e-t-f", 100, -180, 180, 100, 0, 40);
		H_FT_e_t_f.setTitle("electron #theta vs #phi");
		H_FT_e_t_f.setTitleX("#phi (^o)");
		H_FT_e_t_f.setTitleY("#theta (^o)");
		H_FT_e_t_f.setTitle("");
		
		H_FT_e_p_the = new H2F("H-FT-e-p-the", "H-FT-e-p-the", 100, 0, 40, 100, 0, 7);
		H_FT_e_p_the.setTitle("electron p vs #theta (^o)");
		H_FT_e_p_the.setTitleX("#theta (^o)");
		H_FT_e_p_the.setTitleY("p (GeV)");
		H_FT_e_p_the.setTitle("");

		H_FT_e_p_f = new H2F("H-FT-e-p-f", "H-FT-e-p-f", 100, -180, 180, 100, 0, 7);
		H_FT_e_p_f.setTitle("electron p vs #phi");
		H_FT_e_p_f.setTitleX("#phi (^o)");
		H_FT_e_p_f.setTitleY("p (GeV)");
		H_FT_e_p_f.setTitle("");

		H_virphoton = new H1F("H-virphoton", "H-virphoton", 100, 0, 12);
		H_virphoton.setFillColor(LIGHTGREEN);
		H_virphoton.setTitleX("E_#gamma (GeV)");
		H_virphoton.setTitle("");

		hi_rec_fde_theta_vz = new H2F("hi-rec-fde-theta-vz", "hi-rec-fde-theta-vz", 40, -40, 40, 100, 0, 100);
		hi_rec_fde_theta_vz.setTitleX("v_z (cm)");
		hi_rec_fde_theta_vz.setTitleY("#theta (^o)");
		hi_rec_fde_vz = new H1F("hi-rec-fde-vz", "hi-rec-fde-vz", 40, -40, 40);
		hi_rec_fde_vz.setTitleX("v_z (cm)");
		hi_rec_fde_vz.setFillColor(LIGHTGREEN);
		hi_rec_fde_vz.setTitle("");


		// rec electron
		dg_rec_electron = new DataGroup(3, 2);
		dg_rec_electron.addDataSet(H_FT_e_p_the, 0);
		dg_rec_electron.addDataSet(H_FT_e_p_f, 1);
		dg_rec_electron.addDataSet(H_FT_e_t_f, 2);
		//dg_rec_electron.addDataSet(H_virphoton, 3);
		dg_rec_electron.addDataSet(hi_rec_fde_vz, 3);
		dg_rec_electron.addDataSet(hi_rec_fde_theta_vz, 4);

		
		H_FT_W_Q2 = new H2F("H-FT-W-Q2", "H-FT-W-Q2", 100, 0, 5, 100, 0, 12);
		H_FT_W_Q2.setTitle("FT Q^2 vs W");
		H_FT_W_Q2.setTitleX("W ( GeV)");
		H_FT_W_Q2.setTitleY("Q^2 (GeV^2)");
		H_FT_W_Q2.setTitle("");
		
		H_FT_W = new H1F("H-FT-W", "H-FT-W", 100, 0, 5);
		H_FT_W.setTitle("electron W");
		H_FT_W.setTitleX("W (GeV)");
		H_FT_W.setTitleY("count");
		H_FT_W.setFillColor(LIGHTGREEN);
		H_FT_W.setTitle("");


		H_FT_Q2 = new H1F("H-FT-Q2", "H-FT-Q2", 100, 0, 12);
		H_FT_Q2.setFillColor(LIGHTGREEN);
		H_FT_Q2.setTitleY("count");
		H_FT_Q2.setTitleX("Q^2 (GeV^2)");
		H_FT_Q2.setTitle("");



		dg_kinematics = new DataGroup(2,2);
		dg_kinematics.addDataSet(H_FT_Q2, 0);
		dg_kinematics.addDataSet(H_FT_W, 1);
		dg_kinematics.addDataSet(H_FT_W_Q2, 2);
		dg_kinematics.addDataSet(H_virphoton, 3);


	// reconstructed
		// reconstructed first k+ (fast)
		hi_rec_kp1_p_the = new H2F("hi_rec_kp1_p_the", "hi_rec_kp1_p_the", 100, 0, 100, 100, 0, 8);
		hi_rec_kp1_p_the.setTitleX("#theta (^o)");
		hi_rec_kp1_p_the.setTitleY("p (GeV)");
		hi_rec_kp1_p_the.setTitle("");
		hi_rec_kp1_p_phi = new H2F("hi_rec_kp1_p_phi", "hi_rec_kp1_p_phi", 100, -180, 180, 100, 0, 8);
		hi_rec_kp1_p_phi.setTitleX("#phi (^o)");
		hi_rec_kp1_p_phi.setTitleY("p (GeV)");
		hi_rec_kp1_p_phi.setTitle("");
		hi_rec_kp1_the_phi = new H2F("hi_rec_kp1_the_phi", "hi_rec_kp1_the_phi", 100, -180, 180, 100, 0, 100);
		hi_rec_kp1_the_phi.setTitleX("#phi (^o)");
		hi_rec_kp1_the_phi.setTitleY("#theta (^o)");
		hi_rec_kp1_the_phi.setTitle("");
		hi_rec_kp1_theta_vz = new H2F("hi-rec-kp1-theta-vz", "hi-rec-kp1-theta-vz", 40, -40, 40, 100, 0, 100);
		hi_rec_kp1_theta_vz.setTitleX("v_z (cm)");
		hi_rec_kp1_theta_vz.setTitleY("#theta (^o)");
		hi_rec_kp1_theta_vz.setTitle("");
		hi_rec_kp1_vz = new H1F("hi-rec-kp1-vz", "hi-rec-kp1-vz", 40, -40, 40);
		hi_rec_kp1_vz.setTitleX("v_z (cm)");
		hi_rec_kp1_vz.setFillColor(LIGHTGREEN);
		hi_rec_kp1_vz.setTitle("");



		// rec fast kp (kp1)
		dg_rec_kp1 = new DataGroup(3,2);
		dg_rec_kp1.addDataSet(hi_rec_kp1_p_the, 0);
		dg_rec_kp1.addDataSet(hi_rec_kp1_p_phi, 1);
		dg_rec_kp1.addDataSet(hi_rec_kp1_the_phi, 2);
		dg_rec_kp1.addDataSet(hi_rec_kp1_vz, 3);
		dg_rec_kp1.addDataSet(hi_rec_kp1_theta_vz, 4);


		// reconstructed second k+ (slow)
		hi_rec_kp2_p_the = new H2F("hi_rec_kp2_p_the", "hi_rec_kp2_p_the", 100, 0, 100, 100, 0, 8);
		hi_rec_kp2_p_the.setTitleX("#theta (^o)");
		hi_rec_kp2_p_the.setTitleY("p (GeV)");
		hi_rec_kp2_p_the.setTitle("");
		hi_rec_kp2_p_phi = new H2F("hi_rec_kp2_p_phi", "hi_rec_kp2_p_phi", 100, -180, 180, 100, 0, 8);
		hi_rec_kp2_p_phi.setTitleX("#phi (^o)");
		hi_rec_kp2_p_phi.setTitleY("p (GeV)");
		hi_rec_kp2_p_phi.setTitle("");
		hi_rec_kp2_the_phi = new H2F("hi_rec_kp2_the_phi", "hi_rec_kp2_the_phi", 100, -180, 180, 100, 0, 100);
		hi_rec_kp2_the_phi.setTitleX("#phi (^o)");
		hi_rec_kp2_the_phi.setTitleY("#theta (^o)");
		hi_rec_kp2_the_phi.setTitle("");
		hi_rec_kp2_theta_vz = new H2F("hi-rec-kp2-theta-vz", "hi-rec-kp2-theta-vz", 40, -40, 40, 100, 0, 100);
		hi_rec_kp2_theta_vz.setTitleX("v_z (cm)");
		hi_rec_kp2_theta_vz.setTitleY("#theta (^o)");
		hi_rec_kp2_theta_vz.setTitle("");
		hi_rec_kp2_vz = new H1F("hi-rec-kp2-vz", "hi-rec-kp2-vz", 40, -40, 40);
		hi_rec_kp2_vz.setTitleX("v_z (cm)");
		hi_rec_kp2_vz.setFillColor(LIGHTGREEN);
		hi_rec_kp2_vz.setTitle("");


		// rec slow kp (kp2)
		dg_rec_kp2 = new DataGroup(3,2);
		dg_rec_kp2.addDataSet(hi_rec_kp2_p_the, 0);
		dg_rec_kp2.addDataSet(hi_rec_kp2_p_phi, 1);
		dg_rec_kp2.addDataSet(hi_rec_kp2_the_phi, 2);
		dg_rec_kp2.addDataSet(hi_rec_kp2_vz, 3);
		dg_rec_kp2.addDataSet(hi_rec_kp2_theta_vz, 4);


		// two kps 
		hi_rec_kps_dtheta = new H1F("hi-rec-kps-dtheta", "hi-rec-kps-dtheta", 100, -180, 180);
		hi_rec_kps_dtheta.setFillColor(LIGHTGREEN);
		hi_rec_kps_dtheta_slowkp_theta = new H2F("hi-rec-kps-dtheta-slowkp-theta", "hi-rec-kps-dtheta-slowkp-theta", 100, 0, 180, 100, -180, 180);
		hi_rec_kps_dtheta_fastkp_theta = new H2F("hi-rec-kps-dtheta-fastkp-theta", "hi-rec-kps-dtheta-fastkp-theta", 100, 0, 45, 100, -30, 30);
		hi_rec_kps_dtheta_slowkp_phi = new H2F("hi-rec-kps-dtheta-slowkp-phi", "hi-rec-kps-dtheta-slowkp-phi", 100, -180, 180, 100, -30, 30);
		hi_rec_kps_dtheta_fastkp_phi = new H2F("hi-rec-kps-dtheta-fastkp-phi", "hi-rec-kps-dtheta-fastkp-phi", 100, -180, 180, 100, -30, 30);
		hi_rec_kps_dtheta_slowkp_p = new H2F("hi-rec-kps-dtheta-slowkp-p", "hi-rec-kps-dtheta-slowkp-p", 100, 0, 8, 100, -30, 30);
		hi_rec_kps_dtheta_fastkp_p = new H2F("hi-rec-kps-dtheta-fastkp-p", "hi-rec-kps-dtheta-fastkp-p", 100, 0, 8, 100, -30, 30);
		

		dg_rec_kps = new DataGroup(3,2);
		dg_rec_kps.addDataSet(hi_rec_kps_dtheta, 0);
		//dg_rec_kps.addDataSet(hi_rec_kps_dtheta_slowkp_p, 0);
		dg_rec_kps.addDataSet(hi_rec_kps_dtheta_slowkp_theta, 1);
		dg_rec_kps.addDataSet(hi_rec_kps_dtheta_slowkp_phi, 2);
		dg_rec_kps.addDataSet(hi_rec_kps_dtheta_fastkp_p, 3);
		dg_rec_kps.addDataSet(hi_rec_kps_dtheta_fastkp_theta, 4);
		dg_rec_kps.addDataSet(hi_rec_kps_dtheta_fastkp_phi, 5);

		//kps as pips

		hi_efkpfkpfkm_MM_efastkpaspipslowkpfkm = new H2F("hi-efkpfkpfkm-MM-efastkpaspipslowkpfkm", "hi-efkpfkpfkm-MM-efastkpaspipslowkpfkm", 25, 0.9, 1.4, 25, 0.9, 1.4);
		hi_efkpfkpfkm_MM_efastkpaspipslowkpfkm.setTitleX("MM(efK^+(#pi^+)sK^+)");
		hi_efkpfkpfkm_MM_efastkpaspipslowkpfkm.setTitleY("MM(eK^+K^+)");

		hi_efkpfkpfkm_MM_eslowkpaspipfastkpfkm = new H2F("hi-efkpfkpfkm-MM-eslowkpaspipfastkpfkm", "hi-efkpfkpfkm-MM-eslowkpaspipfastkpfkm", 25, 0.9, 1.4, 25, 0.9, 1.4);
		hi_efkpfkpfkm_MM_eslowkpaspipfastkpfkm.setTitleX("MM(esK^+(#pi^+)fK^+)");
		hi_efkpfkpfkm_MM_eslowkpaspipfastkpfkm.setTitleY("MM(eK^+K^+)");

		dg_rec_kpaspip = new DataGroup(3,2);
		dg_rec_kpaspip.addDataSet(hi_efkpfkpfkm_MM_efastkpaspipslowkpfkm, 0);
		dg_rec_kpaspip.addDataSet(hi_efkpfkpfkm_MM_eslowkpaspipfastkpfkm, 1);
		

		// reconstructed  k- 
		hi_rec_km_p_the = new H2F("hi_rec_km_p_the", "hi_rec_km_p_the", 100, 0, 100, 100, 0, 8);
		hi_rec_km_p_the.setTitleX("#theta (^o)");
		hi_rec_km_p_the.setTitleY("p (GeV)");
		hi_rec_km_p_the.setTitle("");
		hi_rec_km_p_phi = new H2F("hi_rec_km_p_phi", "hi_rec_km_p_phi", 100, -180, 180, 100, 0, 8);
		hi_rec_km_p_phi.setTitleX("#phi (^o)");
		hi_rec_km_p_phi.setTitleY("p (GeV)");
		hi_rec_km_p_phi.setTitle("");
		hi_rec_km_the_phi = new H2F("hi_rec_km_the_phi", "hi_rec_km_the_phi", 100, -180, 180, 100, 0, 100);
		hi_rec_km_the_phi.setTitleX("#phi (^o)");
		hi_rec_km_the_phi.setTitleY("#theta (^o)");
		hi_rec_km_the_phi.setTitle("");
		hi_rec_km_theta_vz = new H2F("hi-rec-km-theta-vz", "hi-rec-km-theta-vz", 40, -40, 40, 100, 0, 100);
		hi_rec_km_theta_vz.setTitleX("v_z (cm)");
		hi_rec_km_theta_vz.setTitleY("#theta (^o)");
		hi_rec_km_theta_vz.setTitle("");
		hi_rec_km_vz = new H1F("hi-rec-km-vz", "hi-rec-km-vz", 40, -40, 40);
		hi_rec_km_vz.setTitleX("v_z (cm)");
		hi_rec_km_vz.setFillColor(LIGHTGREEN);
		hi_rec_km_vz.setTitle("");


		// rec km (km);
		dg_rec_km = new DataGroup(3,2);
		dg_rec_km.addDataSet(hi_rec_km_p_the, 0);
		dg_rec_km.addDataSet(hi_rec_km_p_phi, 1);
		dg_rec_km.addDataSet(hi_rec_km_the_phi, 2);
		dg_rec_km.addDataSet(hi_rec_km_vz, 3);
		dg_rec_km.addDataSet(hi_rec_km_theta_vz, 4);


	//counter for no of rec particles 
		hi_pip_counter = new H1F("hi-pip-counter", "hi-pip-counter", 5, 0, 5);
		hi_pip_counter.setFillColor(LIGHTGREEN);
		hi_pim_counter = new H1F("hi-pim-counter", "hi-pim-counter", 5, 0, 5);
		hi_pim_counter.setFillColor(LIGHTGREEN);
		hi_kp_counter = new H1F("hi-kp-counter", "hi-kp-counter", 5, 0, 5);
		hi_kp_counter.setFillColor(LIGHTGREEN);
		hi_km_counter = new H1F("hi-km-counter", "hi-km-counter", 5, 0, 5);
		hi_km_counter.setFillColor(LIGHTGREEN);
		hi_prot_counter = new H1F("hi-prot-counter", "hi-prot-counter", 5, 0, 5);
		hi_prot_counter.setFillColor(LIGHTGREEN);
		hi_fpip_counter = new H1F("hi-fpip-counter", "hi-fpip-counter", 5, 0, 5);
		hi_fpip_counter.setFillColor(LIGHTGREEN);
		hi_fpim_counter = new H1F("hi-fpim-counter", "hi-fpim-counter", 5, 0, 5);
		hi_fpim_counter.setFillColor(LIGHTGREEN);
		hi_fkp_counter = new H1F("hi-fkp-counter", "hi-fkp-counter", 5, 0, 5);
		hi_fkp_counter.setFillColor(LIGHTGREEN);
		hi_fkm_counter = new H1F("hi-fkm-counter", "hi-fkm-counter", 5, 0, 5);
		hi_fkm_counter.setFillColor(LIGHTGREEN);
		hi_fprot_counter = new H1F("hi-fprot-counter", "hi-fprot-counter", 5, 0, 5);
		hi_fprot_counter.setFillColor(LIGHTGREEN);
		hi_cpip_counter = new H1F("hi-cpip-counter", "hi-cpip-counter", 5, 0, 5);
		hi_cpip_counter.setFillColor(LIGHTGREEN);
		hi_cpim_counter = new H1F("hi-cpim-counter", "hi-cpim-counter", 5, 0, 5);
		hi_cpim_counter.setFillColor(LIGHTGREEN);
		hi_ckp_counter = new H1F("hi-ckp-counter", "hi-ckp-counter", 5, 0, 5);
		hi_ckp_counter.setFillColor(LIGHTGREEN);
		hi_ckm_counter = new H1F("hi-ckm-counter", "hi-ckm-counter", 5, 0, 5);
		hi_ckm_counter.setFillColor(LIGHTGREEN);
		hi_cprot_counter = new H1F("hi-cprot-counter", "hi-cprot-counter", 5, 0, 5);
		hi_cprot_counter.setFillColor(LIGHTGREEN);

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

		dg_mm_ekpkpkm = new DataGroup(3, 2);
		dg_mm_ekpkpkm.addDataSet(hi_efkpfkpfkm_MM_efkpfkp, 0);
		dg_mm_ekpkpkm.addDataSet(hi_efkpckpfkm_MM_efkpckp, 1);
		dg_mm_ekpkpkm.addDataSet(hi_eckpckpfkm_MM_eckpckp, 2);
		dg_mm_ekpkpkm.addDataSet(hi_efkpfkpckm_MM_efkpfkp, 3);
		dg_mm_ekpkpkm.addDataSet(hi_efkpckpckm_MM_efkpckp, 4);
		dg_mm_ekpkpkm.addDataSet(hi_eckpckpckm_MM_eckpckp, 5);


		// with lambda events selection in MM(ekpkpkm)
		hi_efkpfkpfkm_MM_efkpfkp_lam_evnt = new H1F("hi-efkpfkpfkm-MM-efkpfkp-lam-evnt", "hi-efkpfkpfkm-MM-efkpfkp-lam-evnt", 30, 1.6, 2.5); //25, 1.6, 2.6//50, 1.6, 2.1 //100, 1.6, 3.1//75, 1.6, 3.1
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

		dg_mm_ekpkpkm_lam_evnt = new DataGroup(3, 2);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_efkpfkpfkm_MM_efkpfkp_lam_evnt, 0);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_efkpckpfkm_MM_efkpckp_lam_evnt, 1);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_eckpckpfkm_MM_eckpckp_lam_evnt, 2);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_efkpfkpckm_MM_efkpfkp_lam_evnt, 3);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_efkpckpckm_MM_efkpckp_lam_evnt, 4);
		dg_mm_ekpkpkm_lam_evnt.addDataSet(hi_eckpckpckm_MM_eckpckp_lam_evnt, 5);



		// IM (kmlambdaconstrained) 
		hi_efkpfkpfkm_IM_kmlambda = new H1F("hi-efkpfkpfkm-IM-kmlambda", "hi-efkpfkpfkm-IM-kmlambda", 30, 1.6, 2.5);//30, 1.6, 2.5 //30, 1.7, 2.6//33, 1.6, 2.6//50, 1.6, 3.1 //50, 1.6, 2.1 //100, 1.6, 3.1//75, 1.6, 3.1
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

		dg_m_kmlambda = new DataGroup(3, 2);
		dg_m_kmlambda.addDataSet(hi_efkpfkpfkm_IM_kmlambda, 0);
		dg_m_kmlambda.addDataSet(hi_efkpckpfkm_IM_kmlambda, 1);
		dg_m_kmlambda.addDataSet(hi_eckpckpfkm_IM_kmlambda, 2);
		dg_m_kmlambda.addDataSet(hi_efkpfkpckm_IM_kmlambda, 3);
		dg_m_kmlambda.addDataSet(hi_efkpckpckm_IM_kmlambda, 4);
		dg_m_kmlambda.addDataSet(hi_eckpckpckm_IM_kmlambda, 5);

		hi_efkpfkpfkm_MM_efkpfkpfkm = new H1F("hi-efkpfkpfkm-MM-efkpfkpfkm", "hi-efkpfkpfkm-MM-efkpfkpfkm", 100, 0.0, 2.4);//25, 0.0, 2.4 //25, 0.9, 1.4
		//hi_efkpfkpfkm_MM_efkpfkpfkm = new H1F("hi-efkpfkpfkm-MM-efkpfkpfkm", "hi-efkpfkpfkm-MM-efkpfkpfkm", 25, 0.9, 1.4);//25, 0.0, 2.4 //25, 0.9, 1.4
		hi_efkpfkpfkm_MM_efkpfkpfkm.setFillColor(LIGHTGREEN);
		hi_efkpckpfkm_MM_efkpckpfkm = new H1F("hi-efkpckpfkm-MM-efkpckpfkm", "hi-efkpckpfkm-MM-efkpckpfkm", 100, 0.0, 2.4);//25, 0.0, 2.4 //25, 0.9, 1.4
		hi_efkpckpfkm_MM_efkpckpfkm.setLineColor(RED);
		hi_efkpckpfkm_MM_efkpckpfkm_nocorn = new H1F("hi-efkpckpfkm-MM-efkpckpfkm-nocorn", "hi-efkpckpfkm-MM-efkpckpfkm-nocorn", 100, 0.0, 2.4);
		hi_efkpckpfkm_MM_efkpckpfkm_nocorn.setFillColor(LIGHTGREEN);
		hi_eckpckpfkm_MM_eckpckpfkm = new H1F("hi-eckpckpfkm-MM-eckpckpfkm", "hi-eckpckpfkm-MM-eckpckpfkm", 100, 0.0, 2.4);
		hi_eckpckpfkm_MM_eckpckpfkm.setFillColor(LIGHTGREEN);
		hi_efkpfkpckm_MM_efkpfkpckm = new H1F("hi-efkpfkpckm-MM-efkpfkpckm", "hi-efkpfkpckm-MM-efkpfkpckm", 100, 0.0, 2.4);
		hi_efkpfkpckm_MM_efkpfkpckm.setFillColor(LIGHTGREEN);
		hi_efkpckpckm_MM_efkpckpckm = new H1F("hi-efkpckpckm-MM-efkpckpckm", "hi-efkpckpckm-MM-efkpckpckm", 100, 0.0, 2.4);
		hi_efkpckpckm_MM_efkpckpckm.setFillColor(LIGHTGREEN);
		hi_eckpckpckm_MM_eckpckpckm = new H1F("hi-eckpckpckm-MM-eckpckpckm", "hi-eckpckpckm-MM-eckpckpckm", 100, 0.0, 2.4);
		hi_eckpckpckm_MM_eckpckp.setFillColor(LIGHTGREEN);

		dg_mm_ekpkpkm1 = new DataGroup(3, 2);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpfkpfkm_MM_efkpfkpfkm, 0);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpckpfkm_MM_efkpckpfkm_nocorn, 1);
		//dg_mm_ekpkpkm1.addDataSet(hi_efkpckpfkm_MM_efkpckpfkm, 1);
		dg_mm_ekpkpkm1.addDataSet(hi_eckpckpfkm_MM_eckpckpfkm, 2);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpfkpckm_MM_efkpfkpckm, 3);
		dg_mm_ekpkpkm1.addDataSet(hi_efkpckpckm_MM_efkpckpckm, 4);
		dg_mm_ekpkpkm1.addDataSet(hi_eckpckpckm_MM_eckpckpckm, 5);

		// scatter plots
		hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm = new H2F("hi-efkpfkpfkm-MM-efkpfkp-MM-efkpfkpfkm", "hi-efkpfkpfkm-MM-efkpfkp-MM-efkpfkpfkm", 100, 0.0, 2.4, 50, 1.0, 3.1); //50, 1.6, 2.1 //100, 1.6, 3.1//75, 1.6, 3.1
		//hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm.setFillColor(LIGHTGREEN);
		hi_efkpckpfkm_MM_efkpckp_MM_efkpckpfkm = new H2F("hi-efkpckpfkm-MM-efkpckp-MM-efkpckpfkm", "hi-efkpckpfkm-MM-efkpckp-MM-efkpckpfkm", 100, 0.0, 2.4, 50, 1.0, 3.1);
		//hi_efkpckpfkm_MM_efkpckp_MM_efkpckpfkm.setFillColor(LIGHTGREEN);
		hi_eckpckpfkm_MM_eckpckp_MM_eckpckpfkm = new H2F("hi-eckpckpfkm-MM-eckpckp-MM-eckpckpfkm", "hi-eckpckpfkm-MM-eckpckp-MM-eckpckpfkm", 100, 0.0, 2.4, 50, 1.0, 3.1);
		//hi_eckpckpfkm_MM_eckpckp_MM_eckpckpfkm.setFillColor(LIGHTGREEN);
		hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm = new H2F("hi-efkpfkpckm-MM-efkpfkp-MM-efkpfkpckm", "hi-efkpfkpckm-MM-efkpfkp-MM-efkpfkpckm", 100, 0.0, 2.4, 50, 1.0, 3.1);
		//hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm.setFillColor(LIGHTGREEN);
		hi_efkpckpckm_MM_efkpckp_MM_efkpckpckm = new H2F("hi-efkpckpckm-MM-efkpckp-MM-efkpckpckm", "hi-efkpckpckm-MM-efkpckp-MM-efkpckpckm", 100, 0.0, 2.4, 50, 1.0, 3.1);
		//hi_efkpckpckm_MM_efkpckp_MM_efkpckpckm.setFillColor(LIGHTGREEN);
		hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm = new H2F("hi-eckpckpckm-MM-eckpckp-MM-eckpckpckm", "hi-eckpckpckm-MM-eckpckp-MM-eckpckpckm", 100, 0.0, 2.4, 50, 1.0, 3.1);
		//hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm.setFillColor(LIGHTGREEN);

		dg_mm_scatter = new DataGroup(3, 2);
		dg_mm_scatter.addDataSet(hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm, 0);
		dg_mm_scatter.addDataSet(hi_efkpckpfkm_MM_efkpckp_MM_efkpckpfkm, 1);
		dg_mm_scatter.addDataSet(hi_eckpckpfkm_MM_eckpckp_MM_eckpckpfkm, 2);
		dg_mm_scatter.addDataSet(hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm, 3);
		dg_mm_scatter.addDataSet(hi_efkpckpckm_MM_efkpckp_MM_efkpckpckm, 4);
		dg_mm_scatter.addDataSet(hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm, 5);
		



	
	//ekpkpreqkm	reconstructed
		
		hi_ekpkpkm_MM_ekpkp = new H1F("hi-ekpkpkm-MM-ekpkp", "hi-ekpkpkm-MM-ekpkp", 50, 1.6, 2.1);
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
    	hi_ekpkpkm_MM_ekpkp_nocorr = new H1F("hi-ekpkpkm-MM-ekpkp-nocorr", "hi-ekpkpkm-MM-ekpkp-nocorr", 50, 1.6, 2.1);
		hi_ekpkpkm_MM_ekpkp_nocorr.setTitle("MM");
		hi_ekpkpkm_MM_ekpkp_nocorr.setTitleX("MM(eK^+K^+) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkp_nocorr.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkp_nocorr.setFillColor(LIGHTGREEN);

		fn_xi_no_corr = new F1D("fn-xi-no-corr", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		fn_xi_no_corr.setParameter(0, 0);
    	fn_xi_no_corr.setParameter(1, 1);
    	fn_xi_no_corr.setParameter(2, 0.2);
    	fn_xi_no_corr.setLineWidth(2);
    	fn_xi_no_corr.setLineColor(2);
    	fn_xi_no_corr.setOptStat("1111");

    //ekpkpkm generated
    	
    	hi_mc_ekpkpkm_mm_ekpkp = new H1F("hi-mc-ekpkpkm-mm-ekpkp", "hi-mc-ekpkpkm-mm-ekpkp", 50, 1.6, 2.1);
		hi_mc_ekpkpkm_mm_ekpkp.setTitle("MM");
		hi_mc_ekpkpkm_mm_ekpkp.setTitleX("MM(eK^+K^+) [GeV/c^2]");
		hi_mc_ekpkpkm_mm_ekpkp.setTitleY("Events/[10 MeV/c^2]");
		hi_mc_ekpkpkm_mm_ekpkp.setFillColor(LIGHTGREEN);

		f1_mc_xi = new F1D("f1-mc-xi", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		//f1_xi.setParameter(0, 0);
    	//f1_xi.setParameter(1, 1);
    	//f1_xi.setParameter(2, 0.2);
    	f1_mc_xi.setLineWidth(2);
    	f1_mc_xi.setLineColor(2);
    	f1_mc_xi.setOptStat("1111");

    	hi_mc_ekpkpkm_mm_ekpkpkm = new H1F("hi-mc-ekpkpkm-mm-ekpkpkm", "hi-mc-ekpkpkm-mm-ekpkpkm", 50, 0.9, 1.4);
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
		hi_ekpkpkm_MM_ekpkpkm = new H1F("hi-ekpkpkm-MM-ekpkpkm", "hi-ekpkpkm-MM-ekpkpkm", 50, 0.9, 1.4);
		hi_ekpkpkm_MM_ekpkpkm.setTitle("MM");
		hi_ekpkpkm_MM_ekpkpkm.setTitleX("MM(eK^+K^+K^-) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm.setFillColor(LIGHTGREEN);

		// for one with common vertex correction for electron
		hi_ekpkpkm_MM_ekpkpkm_nocorr = new H1F("hi-ekpkpkm-MM-ekpkpkm-nocorr", "hi-ekpkpkm-MM-ekpkpkm-nocorr", 50, 0.9, 1.4);
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitle("MM");
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitleX("MM(eK^+K^+K^-) [GeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_MM_ekpkpkm_nocorr.setFillColor(LIGHTGREEN);

		hi_ekpkpkm_IM_kmlambda = new H1F("hi-ekpkpkm-IM-kmlambda", "hi-ekpkpkm-IM-kmlambda", 50, 1.6, 2.1);
		hi_ekpkpkm_IM_kmlambda.setTitle("M");
		hi_ekpkpkm_IM_kmlambda.setTitleX("M(K^-#Lambda) [GeV/c^2]");
		hi_ekpkpkm_IM_kmlambda.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_IM_kmlambda.setFillColor(LIGHTGREEN);

		fn_im_kmlambda = new F1D("fn-im-kmlambda", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		//f1_xi.setParameter(0, 0);
    	//f1_xi.setParameter(1, 1);
    	//f1_xi.setParameter(2, 0.2);
    	fn_im_kmlambda.setLineWidth(2);
    	fn_im_kmlambda.setLineColor(2);
    	fn_im_kmlambda.setOptStat("1111");

		hi_ekpkpkm_IM_kmsigma = new H1F("hi-ekpkpkm-IM-kmsigma", "hi-ekpkpkm-IM-kmsigma", 50, 1.6, 2.1);
		hi_ekpkpkm_IM_kmsigma.setTitle("M");
		hi_ekpkpkm_IM_kmsigma.setTitleX("M(K^-#Sigma) [GeV/c^2]");
		hi_ekpkpkm_IM_kmsigma.setTitleY("Events/[10 MeV/c^2]");
		hi_ekpkpkm_IM_kmsigma.setFillColor(LIGHTGREEN);

		fn_im_kmsigma = new F1D("fn-im-kmsigma", "[amp]*gaus(x,[mean],[sigma])", 1.7, 1.9);
		//f1_xi.setParameter(0, 0);
    	//f1_xi.setParameter(1, 1);
    	//f1_xi.setParameter(2, 0.2);
    	fn_im_kmsigma.setLineWidth(2);
    	fn_im_kmsigma.setLineColor(2);
    	fn_im_kmsigma.setOptStat("1111");


		f1_lambda = new F1D("f1-lambda", "[amp]*gaus(x,[mean],[sigma])", 1.04, 1.14);
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

		dg_mm = new DataGroup(3,3);
		dg_mm.addDataSet(hi_ekpkp_MM_req_pim, 0);
		dg_mm.addDataSet(hi_efkpckp_MM_req_pim, 1);
		dg_mm.addDataSet(hi_efkpckp_MM_req_cpim, 2);
		dg_mm.addDataSet(hi_ekpkpkm_MM_ekpkp, 3);
		dg_mm.addDataSet(hi_efkpckp_MM_req_km, 4);
		dg_mm.addDataSet(hi_efkpckp_MM_req_ckm, 5);
		dg_mm.addDataSet(hi_ekpkpkmprot_MM2, 6);



		////////// PARTICLES vertex/pathlength/TOF/VertexTimeVsP/BetaVsP Plots (saperated according to CD and FD phase space)///////////////
		
		//FD particle vertex
		hi_pim_vz = new H1F("H-pim-vz", "H-pim-vz", 100, -20, 20);
		hi_pim_vz.setFillColor(LIGHTGREEN);
		hi_pip_vz = new H1F("H-pip-vz", "H-pip-vz", 100, -20, 20);
		hi_pip_vz.setFillColor(LIGHTGREEN);
		hi_kp_vz = new H1F("H-kp-vz", "H-kp-vz", 100, -20, 20);//events/0.4 cm
		hi_kp_vz.setFillColor(LIGHTGREEN);
		hi_km_vz = new H1F("H-km-vz", "H-km-vz", 100, -20, 20);
		hi_km_vz.setFillColor(LIGHTGREEN);
		hi_prot_vz = new H1F("H-prot-vz", "H-prot-vz", 100, -20, 20);
		hi_prot_vz.setFillColor(LIGHTGREEN);
		
		//FD particle path
		hi_prot_FTOF1b_path = new H1F("H-prot-FTOF1b-path", "H-prot-FTOF1b-path", 200, 400, 900);
		hi_prot_FTOF1b_path.setFillColor(LIGHTGREEN);
		hi_pip_FTOF1b_path = new H1F("H-pip-FTOF1b-path", "H-pip-FTOF1b-path", 200, 400, 900);
		hi_pip_FTOF1b_path.setFillColor(LIGHTGREEN);
		hi_pim_FTOF1b_path = new H1F("H-pim-FTOF1b-path", "H-pim-FTOF1b-path", 200, 400, 900);
		hi_pim_FTOF1b_path.setFillColor(LIGHTGREEN);
		hi_kp_FTOF1b_path = new H1F("H-kp-FTOF1b-path", "H-kp-FTOF1b-path", 200, 400, 900);
		hi_kp_FTOF1b_path.setFillColor(LIGHTGREEN);
		hi_km_FTOF1b_path = new H1F("H-km-FTOF1b-path", "H-km-FTOF1b-path", 200, 400, 900);
		hi_km_FTOF1b_path.setFillColor(LIGHTGREEN);
		//FD particle time
		hi_prot_FTOF1b_t = new H1F("H-prot-FTOF1b-t", "H-prot-FTOF1b-t", 200, 0, 80);
		hi_prot_FTOF1b_t.setFillColor(LIGHTGREEN);
		hi_pip_FTOF1b_t = new H1F("H-pip-FTOF1b-t", "H-pip-FTOF1b-t", 200, 0, 80);
		hi_pip_FTOF1b_t.setFillColor(LIGHTGREEN);
		hi_pim_FTOF1b_t = new H1F("H-pim-FTOF1b-t", "H-pim-FTOF1b-t", 200, 0, 80);
		hi_pim_FTOF1b_t.setFillColor(LIGHTGREEN);
		hi_kp_FTOF1b_t = new H1F("H-kp-FTOF1b-t", "H-kp-FTOF1b-t", 200, 0, 80);
		hi_kp_FTOF1b_t.setFillColor(LIGHTGREEN);
		hi_km_FTOF1b_t = new H1F("H-km-FTOF1b-t", "H-km-FTOF1b-t", 200, 0, 80);
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
		hi_pimc_vz = new H1F("H-pimc-vz", "H-pimc-vz", 100, -20, 20);
		hi_pimc_vz.setFillColor(LIGHTGREEN);
		hi_pipc_vz = new H1F("H-pipc-vz", "H-pipc-vz", 100, -20, 20);
		hi_pipc_vz.setFillColor(LIGHTGREEN);
		hi_kpc_vz = new H1F("H-kpc-vz", "H-kpc-vz", 100, -20, 20);
		hi_kpc_vz.setFillColor(LIGHTGREEN);
		hi_kmc_vz = new H1F("H-km-vz", "H-kmc-vz", 100, -20, 20);
		hi_kmc_vz.setFillColor(LIGHTGREEN);
		hi_protc_vz = new H1F("H-protc-vz", "H-protc-vz", 100, -20, 20);
		hi_protc_vz.setFillColor(LIGHTGREEN);
		
		//CD particle path
		hi_prot_CTOF_path = new H1F("H-prot-CTOF-path", "H-prot-CTOF-path", 200, 0, 100);
		hi_prot_CTOF_path.setFillColor(LIGHTGREEN);
		hi_pip_CTOF_path = new H1F("H-pip-CTOF-path", "H-pip-CTOF-path", 200, 0, 100);
		hi_pip_CTOF_path.setFillColor(LIGHTGREEN);
		hi_pim_CTOF_path = new H1F("H-pim-CTOF-path", "H-pim-CTOF-path", 200, 0, 100);
		hi_pim_CTOF_path.setFillColor(LIGHTGREEN);
		hi_kp_CTOF_path = new H1F("H-kp-CTOF-path", "H-kp-CTOF-path", 200, 0, 100);
		hi_kp_CTOF_path.setFillColor(LIGHTGREEN);
		hi_km_CTOF_path = new H1F("H-km-CTOF-path", "H-km-CTOF-path", 200, 0, 100);
		hi_km_CTOF_path.setFillColor(LIGHTGREEN);
		
		//CD particle time
		hi_prot_CTOF_t = new H1F("H-prot-CTOF-t", "H-prot-CTOF-t", 200, -5, 20);
		hi_prot_CTOF_t.setFillColor(LIGHTGREEN);
		hi_pip_CTOF_t = new H1F("H-pip-CTOF-t", "H-pip-CTOF-t", 200, -5, 20);
		hi_pip_CTOF_t.setFillColor(LIGHTGREEN);
		hi_pim_CTOF_t = new H1F("H-pim-CTOF-t", "H-pim-CTOF-t", 200, -5, 20);
		hi_pim_CTOF_t.setFillColor(LIGHTGREEN);
		hi_kp_CTOF_t = new H1F("H-kp-CTOF-t", "H-kp-CTOF-t", 200, -5, 20);
		hi_kp_CTOF_t.setFillColor(LIGHTGREEN);
		hi_km_CTOF_t = new H1F("H-km-CTOF-t", "H-km-CTOF-t", 200, -5, 20);
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
		//hi_FD_pos_beta_mom = new H2F("H-FD-pos-beta-mom", "H-FD-pos-beta-mom", 80, 0, 8, 50, 0.3, 1.1);
		hi_FD_pos_beta_mom = new H2F("H-FD-pos-beta-mom", "H-FD-pos-beta-mom", 250, 0.5, 5, 200, 0.85, 1.05); //80, 0.4, 5, 50, 0.85, 1.05
		hi_FD_pos_beta_mom.setTitle("POS  #beta vs mom");
		hi_FD_pos_beta_mom.setTitleX("p (GeV)");
		hi_FD_pos_beta_mom.setTitleY("#beta");
		//hi_FD_neg_beta_mom = new H2F("H-FD-neg-beta-mom", "H-FD-neg-beta-mom", 80, 0, 8, 50, 0.3, 1.1);
		hi_FD_neg_beta_mom = new H2F("H-FD-neg-beta-mom", "H-FD-neg-beta-mom", 250, 0.5, 5, 200, 0.85, 1.05); //80, 0.4, 5, 50, 0.85, 1.05
		hi_FD_neg_beta_mom.setTitle("NEG  #beta vs mom");
		hi_FD_neg_beta_mom.setTitleX("p (GeV)");
		hi_FD_neg_beta_mom.setTitleY("#beta");
		hi_FD_neutral_beta_mom = new H2F("H-FD-neutral-beta-mom", "H-FD-neutral-beta-mom", 80, 0, 8, 50, 0.3, 1.1);
		hi_FD_neutral_beta_mom.setTitle("NEUTRAL  #beta vs mom");
		hi_FD_neutral_beta_mom.setTitleX("p (GeV)");
		hi_FD_neutral_beta_mom.setTitleY("FTB #beta");
		hi_fd_pos_mass = new H1F("hi-fd-pos-mass", "hi-fd-pos-mass", 150, -0.5, 4.5);
		hi_fd_pos_mass.setFillColor(LIGHTGREEN);
		hi_FD_pos_mass_mom = new H2F("H-FD-pos-mass-mom", "H-FD-pos-mass-mom", 100, 0, 7, 150, -0.5, 4.5);
		hi_FD_pos_mass_mom.setTitle("POS Mass^2 vs mom");
		hi_FD_pos_mass_mom.setTitleX("p (GeV)");
		hi_FD_pos_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_fd_neg_mass = new H1F("hi-fd-neg-mass", "hi-fd-neg-mass", 150, -0.5, 2.0);
		hi_fd_neg_mass.setFillColor(LIGHTGREEN);
		hi_FD_neg_mass_mom = new H2F("H-FD-neg-mass-mom", "H-FD-neg-mass-mom", 100, 0, 7, 150, -0.5, 2.0);
		hi_FD_neg_mass_mom.setTitle("NEG Mass^2 vs mom");
		hi_FD_neg_mass_mom.setTitleX("p (GeV)");
		hi_FD_neg_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_FD_neutral_mass_mom = new H2F("H-FD-neutral-mass-mom", "H-FD-neutral-mass-mom", 100, 0, 7, 150, -0.5, 2.0);
		hi_FD_neutral_mass_mom.setTitle("NEUTRAL Mass^2 vs mom");
		hi_FD_neutral_mass_mom.setTitleX("p (GeV)");
		hi_FD_neutral_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_FD_pos_mass_the = new H2F("H-FD-pos-mass-the", "H-FD-pos-mass-the", 100, 0, 80, 100, -0.5, 4.5);
		hi_FD_pos_mass_the.setTitle("POS Mass^2 vs #theta");
		hi_FD_pos_mass_the.setTitleX("#theta (^o)");
		hi_FD_pos_mass_the.setTitleY("M^2 (GeV^2");
		hi_FD_neg_mass_the = new H2F("H-FD-neg-mass-the", "H-FD-neg-mass-the", 100, 0, 80, 100, -0.5, 2);
		hi_FD_neg_mass_the.setTitle("NEG Mass^2 vs #theta");
		hi_FD_neg_mass_the.setTitleX("#theta (^o)");
		hi_FD_neg_mass_the.setTitleY("M^2 (GeV^2");
		hi_FD_neutral_mass_the = new H2F("H-FD-neutral-mass-the", "H-FD-neutral-mass-the", 100, 0, 80, 100, -0.5, 2);
		hi_FD_neutral_mass_the.setTitle("NEUTRAL Mass^2 vs #theta");
		hi_FD_neutral_mass_the.setTitleX("#theta (^o)");
		hi_FD_neutral_mass_the.setTitleY("M^2 (GeV^2");
		hi_FD_pos_mass_phi = new H2F("H-FD-pos-mass-phi", "H-FD-pos-mass-phi", 100, -180, 180, 100, -0.5, 4.5);
		hi_FD_pos_mass_phi.setTitle("POS Mass^2 vs #phi");
		hi_FD_pos_mass_phi.setTitleX("#phi (^o)");
		hi_FD_pos_mass_phi.setTitleY("M^2 (GeV^2");
		hi_FD_neg_mass_phi = new H2F("H-FD-neg-mass-phi", "H-FD-neg-mass-phi", 100, -180, 180, 100, -0.5, 2);
		hi_FD_neg_mass_phi.setTitle("NEG Mass^2 vs #phi");
		hi_FD_neg_mass_phi.setTitleX("#phi (^o)");
		hi_FD_neg_mass_phi.setTitleY("M^2 (GeV^2");
		hi_FD_neutral_mass_phi = new H2F("H-FD-neutral-mass-phi", "H-FD-neutral-mass-phi", 100, -180, 180, 100, -0.5, 2);
		hi_FD_neutral_mass_phi.setTitle("NEUTRAL Mass^2 vs #phi");
		hi_FD_neutral_mass_phi.setTitleX("#phi (^o)");
		hi_FD_neutral_mass_phi.setTitleY("M^2 (GeV^2");
		
		
		// CD particle beta vs mom by charge
		hi_CD_pos_beta_mom = new H2F("H-CD-pos-beta-mom", "H-CD-pos-beta-mom", 50, 0.1, 2.0, 50, 0.1, 1.1);
		hi_CD_pos_beta_mom.setTitle("POS  #beta vs mom");
		hi_CD_pos_beta_mom.setTitleX("p (GeV)");
		hi_CD_pos_beta_mom.setTitleY("FTB #beta");
		hi_CD_neg_beta_mom = new H2F("H-CD-neg-beta-mom", "H-CD-neg-beta-mom", 50, 0.1, 2.0, 50, 0.1, 1.1);
		hi_CD_neg_beta_mom.setTitle("NEG  #beta vs mom");
		hi_CD_neg_beta_mom.setTitleX("p (GeV)");
		hi_CD_neg_beta_mom.setTitleY("FTB #beta");
		hi_CD_neutral_beta_mom = new H2F("H-CD-neutral-beta-mom", "H-CD-neutral-beta-mom", 50, 0, 2.0, 50, 0.1, 1.1);
		hi_CD_neutral_beta_mom.setTitle("NEUTRAL  #beta vs mom");
		hi_CD_neutral_beta_mom.setTitleX("p (GeV)");
		hi_CD_neutral_beta_mom.setTitleY("FTB #beta");
		hi_cd_pos_mass = new H1F("hi-cd-pos-mass", "hi-cd-pos-mass", 150, -0.5, 4.5);
		hi_cd_pos_mass.setFillColor(LIGHTGREEN);
		hi_CD_pos_mass_mom = new H2F("H-CD-pos-mass-mom", "H-CD-pos-mass-mom", 50, 0, 2.0, 150, -0.5, 4.5);
		hi_CD_pos_mass_mom.setTitle("POS Mass^2 vs mom");
		hi_CD_pos_mass_mom.setTitleX("p (GeV)");
		hi_CD_pos_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_cd_neg_mass = new H1F("hi-fd-neg-mass", "hi-fd-neg-mass", 150, -0.5, 2.0);
		hi_cd_neg_mass.setFillColor(LIGHTGREEN);
		hi_CD_neg_mass_mom = new H2F("H-CD-neg-mass-mom", "H-CD-neg-mass-mom", 50, 0, 2.0, 150, -0.5, 2.0);
		hi_CD_neg_mass_mom.setTitle("NEG Mass^2 vs mom");
		hi_CD_neg_mass_mom.setTitleX("p (GeV)");
		hi_CD_neg_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_CD_neutral_mass_mom = new H2F("H-CD-neutral-mass-mom", "H-CD-neutral-mass-mom", 50, 0, 2.0, 150, -0.5, 2.0);
		hi_CD_neutral_mass_mom.setTitle("NEUTRAL Mass^2 vs mom");
		hi_CD_neutral_mass_mom.setTitleX("p (GeV)");
		hi_CD_neutral_mass_mom.setTitleY("M^2 (GeV^2)");
		hi_CD_pos_mass_the = new H2F("H-CD-pos-mass-the", "H-CD-pos-mass-the", 100, 30, 130, 100, -0.5, 4.5);
		hi_CD_pos_mass_the.setTitle("POS Mass^2 vs #theta");
		hi_CD_pos_mass_the.setTitleX("#theta (^o)");
		hi_CD_pos_mass_the.setTitleY("M^2 (GeV^2");
		hi_CD_neg_mass_the = new H2F("H-CD-neg-mass-the", "H-CD-neg-mass-the", 100, 30, 130, 100, -0.5, 2);
		hi_CD_neg_mass_the.setTitle("NEG Mass^2 vs #theta");
		hi_CD_neg_mass_the.setTitleX("#theta (^o)");
		hi_CD_neg_mass_the.setTitleY("M^2 (GeV^2");
		hi_CD_neutral_mass_the = new H2F("H-CD-neutral-mass-the", "H-CD-neutral-mass-the", 100, 30, 130, 100, -0.5, 2);
		hi_CD_neutral_mass_the.setTitle("NEUTRAL Mass^2 vs #theta");
		hi_CD_neutral_mass_the.setTitleX("#theta (^o)");
		hi_CD_neutral_mass_the.setTitleY("M^2 (GeV^2");
		hi_CD_pos_mass_phi = new H2F("H-CD-pos-mass-phi", "H-CD-pos-mass-phi", 150, -180, 180, 100, -0.5, 4.5);
		hi_CD_pos_mass_phi.setTitle("POS Mass^2 vs #phi");
		hi_CD_pos_mass_phi.setTitleX("#phi (^o)");
		hi_CD_pos_mass_phi.setTitleY("M^2 (GeV^2");
		hi_CD_neg_mass_phi = new H2F("H-CD-neg-mass-phi", "H-CD-neg-mass-phi", 150, -180, 180, 100, -0.5, 2);
		hi_CD_neg_mass_phi.setTitle("NEG Mass^2 vs #phi");
		hi_CD_neg_mass_phi.setTitleX("#phi (^o)");
		hi_CD_neg_mass_phi.setTitleY("M^2 (GeV^2");
		hi_CD_neutral_mass_phi = new H2F("H-CD-neutral-mass-phi", "H-CD-neutral-mass-phi", 150, -180, 180, 100, -0.5, 2);
		hi_CD_neutral_mass_phi.setTitle("NEUTRAL Mass^2 vs #phi");
		hi_CD_neutral_mass_phi.setTitleX("#phi (^o)");
		hi_CD_neutral_mass_phi.setTitleY("M^2 (GeV^2");

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

		//CD protons and pips
		hi_cd_prot_p = new H1F("hi_cd_prot_p", "hi-cd-prot-p", 50, 0, 2.0);
		hi_cd_prot_p.setTitleX("p (GeV)");
		hi_cd_prot_p.setTitleY("count");
		hi_cd_prot_p.setFillColor(LIGHTGREEN);
		hi_cd_pip_p = new H1F("hi_cd_prot_p", "hi-cd-pip-p", 50, 0, 2.0);
		hi_cd_pip_p.setTitleX("p (GeV)");
		hi_cd_pip_p.setTitleY("count");
		hi_cd_pip_p.setFillColor(LIGHTGREEN);
		hi_cd_prot_theta = new H1F("hi_cd_prot_theta", "hi-cd-prot-theta", 150, 30, 130);
		hi_cd_prot_theta.setTitleX("#theta (^o)");
		hi_cd_prot_theta.setTitleY("count");
		hi_cd_prot_theta.setFillColor(LIGHTGREEN);
		hi_cd_pip_theta = new H1F("hi_cd_prot_theta", "hi-cd-pip-theta", 150, 30, 130);
		hi_cd_pip_theta.setTitleX("#theta (^o)");
		hi_cd_pip_theta.setTitleY("count");
		hi_cd_pip_theta.setFillColor(LIGHTGREEN);
		hi_cd_prot_p_theta = new H2F("hi-cd-prot-p-theta", "hi-cd-prot-p-theta", 150, 30, 130, 50, 0, 2.0);
		hi_cd_prot_p_theta.setTitleX("#theta (^o)");
		hi_cd_prot_p_theta.setTitleY("p (GeV)");
		hi_cd_pip_p_theta = new H2F("hi-cd-pip-p-theta", "hi-cd-pip-p-theta", 150, 30, 130, 50, 0, 2.0);
		hi_cd_pip_p_theta.setTitleX("#theta (^o)");
		hi_cd_pip_p_theta.setTitleY("p (GeV)");

		hi_W_cd_pro_the = new H2F("hi-W-cd-pro-the", "hi-W-cd-pro-the", 100, 0, 5, 100, 30, 130);
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


		
	} // end of ft_ana()


///*
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




		/*
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

	    	//dg_rec_y.getH1F("hi_mc_ekpkpkm_mm_ekpkp").fill(xi.mass());
	    	//dg_rec_y.getH1F("hi_mc_ekpkpkm_mm_ekpkpkm").fill(lambda.mass()); 


		}*/


 	}

	//*/

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

	//public int makeFDElectron(DataBank bank){
	public int makeFDElectron(DataBank bank, HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank, HipoDataBank cc_Bank) {
		int NFDElec = 0;
		fdels = new ArrayList<Particle>();
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
			double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));

			if (pid == 11 && inDC && bank.getShort("status", k) < 0.0 && q == -1 && e_mom >= 1.0){//

				//to check
				hi_pim_vz.fill(vz);
			}
			

			if (pid == 11 && inDC && bank.getShort("status", k) < 0.0 && q == -1 && e_mom >= 1.0 && vz >= -10.0 && vz <= 2.0 ){ //&& e_the >= 5.0 && e_the <= 35.0
				NFDElec++;
				Particle Vfdel = new Particle(pid, px, py, pz, vx, vy, vz);

				fdels.add(Vfdel);
			}
			//Timothy's analysis fitter
            //GenericKinematicFitter research_fitter = new analysis_fitter(Eb);

            if (NFDElec == 1 ){ //sector 1 hit
			//if (NFDElec == 1 
			//	&& research_fitter.particle_test(k, bank) 
			//	&& research_fitter.electron_test(k, p, vz, bank, cal_Bank, track_Bank, traj_Bank, run_Bank, cc_Bank)) {

				found_eFD = true;
				double eft = Math.sqrt(e_mom * e_mom + PDGDatabase.getParticleById(11).mass()*PDGDatabase.getParticleById(11).mass());
				e_phi = (float) Math.toDegrees(Math.atan2(py, px));
				Vfde = new Particle(pid, px, py, pz, vx, vy, vz);
				//Ve = new LorentzVector(px, py, pz, eft);

				if (runType == "mc"){
					//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
					//Random rand = new Random();
    				//double smearFactor = rand.nextGaussian()/100;
    				double smearedP = Vfde.p()*(1+smearFactor);
					Vfde.setP(smearedP);

					//Giovanni smearing function
					//LorentzVector lv_smeared = new LorentzVector();
					//lv_smeared = smear(Vfde.vector());
					//Vfde.setVector(lv_smeared, Vfde.vertex());

				}

				Ve = Vfde.vector();

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

	public void makeOthers(DataBank recbank, DataBank recSCBank, HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank, HipoDataBank cc_Bank) {
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

		// momentum smearing factor for MC reconstruction
		//Random rand = new Random();
    	//double smearFactor = rand.nextGaussian()/100;
		
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
							//pip_FTOF1b_vt = pip_FTOF1b_t - pip_FTOF1b_path / (pip_beta * 29.98f) - STT - pip_vz/ (pip_beta * 29.98f);
							//pip_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							float pip_FTOF1b_TOF = (float) pip_FTOF1b_t - STT - pip_vz/ (pip_beta * 29.98f);
							//difference in the measured and computed vertex times 
							pip_FTOF1b_vt = (float) pip_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass())));
							if (Math.abs(pip_FTOF1b_vt) < 0.5 && beta > 0.4 && beta < 1.05 && mom > 0.4 && mom < Eb) { //0.4
								nfdwithcutpip++;
								//Vpip = new LorentzVector(pip_px*pipCor, pip_py*pipCor, pip_pz*pipCor, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								//Vpip = new LorentzVector(pip_px, pip_py, pip_pz, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								Vpip = new Particle(pid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);
								fpips.add(Vpip);
								pips.add(Vpip);
								//hi_pip_vt_p.fill(pip_FTOF1b_vt, pip_mom);
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
							float pip_CTOF_TOF = (float) pip_CTOF_t - STT - pip_vz/ (pipc_beta * 29.98f);
							//difference in the measured and computed vertex times 
							pip_CTOF_vt = (float) pip_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass())));
							if (Math.abs(pip_CTOF_vt) < 0.4 && beta > 0.2 && beta < 1.05 && mom > 0.2 && mom < 3.0) {
								ncdwithcutpip++;
								//Vpipc = new LorentzVector(pip_px*pipCor, pip_py*pipCor, pip_pz*pipCor, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
								Vpipc = new Particle(pid, pip_px, pip_py, pip_pz, pip_vx, pip_vy, pip_vz);	
								pips.add(Vpipc);
								cpips.add(Vpipc);
								//hi_pipc_vt_p.fill(pip_FTOF1b_vt, pip_mom);
							}
							hi_pipc_vt_p.fill(pip_mom, pip_CTOF_vt);
							hi_pipc_vz.fill(pip_vz);
							hi_pip_CTOF_t.fill(pip_CTOF_t - STT - pip_vz/ (pipc_beta * 29.98f));
							hi_pip_CTOF_path.fill(pip_CTOF_path);
							
							
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
					if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut ) {	//&& research_fitter.dc_fiducial_cut(k, recbank, track_Bank, traj_Bank, run_Bank)
						if (recSCBank.getShort("pindex", r) == pim_part_ind && recSCBank.getByte("layer", r) == 2) {
							pim_FTOF_pad1b = recSCBank.getShort("component", r);
							pim_FTOF1b_t = recSCBank.getFloat("time", r);
							pim_FTOF1b_path = recSCBank.getFloat("path", r);
							float pim_beta = pim_mom / (float) Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass() * PDGDatabase.getParticleById(-211).mass());
							//pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_ftb_beta * 29.98f) - STT - pim_vz/ (pim_ftb_beta * 29.98f);
							//pim_FTOF1b_vt = pim_FTOF1b_t - pim_FTOF1b_path / (pim_beta * 29.98f) - STT - pim_vz/ (pim_beta * 29.98f);
							//pim_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							float pim_FTOF1b_TOF = (float) pim_FTOF1b_t - STT - pim_vz/ (pim_beta * 29.98f);
							//difference in the measured and computed vertex times 
							pim_FTOF1b_vt = (float) pim_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass())));
							
							if (Math.abs(pim_FTOF1b_vt) < 0.5 && beta > 0.4 && beta < 1.05 && mom > 0.4 && mom < Eb) {
								//Vpim = new LorentzVector(pim_px*pimCor, pim_py*pimCor, pim_pz*pimCor, Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass()));
								Vpim = new Particle(pid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);
								pims.add(Vpim);
								fpims.add(Vpim);
								//hi_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);
							}
							hi_pim_vt_p.fill(pim_mom, pim_FTOF1b_vt);
							//hi_pim_vz.fill(pim_vz);
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
							//pim_CTOF_vt = pim_CTOF_t - pim_CTOF_path / (pimc_beta * 29.98f) - STT - pim_vz/ (pimc_beta * 29.98f);
							//pim_CTOF_vt = STT - recFTbank.getFloat("vt", k);
							float pim_CTOF_TOF = (float) pim_CTOF_t - STT - pim_vz/ (pimc_beta * 29.98f);
							//difference in the measured and computed vertex times 
							pim_CTOF_vt = (float) pim_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass())));
							if (Math.abs(pim_CTOF_vt) < 0.4 && beta > 0.2 && beta < 1.05 && mom > 0.2 && mom < 3.0) {
								//Vpimc = new LorentzVector(pim_px*pimCor, pim_py*pimCor, pim_pz*pimCor, Math.sqrt(pim_mom * pim_mom + PDGDatabase.getParticleById(-211).mass()*PDGDatabase.getParticleById(-211).mass()));
								Vpimc = new Particle(pid, pim_px, pim_py, pim_pz, pim_vx, pim_vy, pim_vz);		
								pims.add(Vpimc);
								cpims.add(Vpimc);
								//hi_pim_vt_p.fill(pim_FTOF1b_vt, pim_mom);
							}
							hi_pimc_vt_p.fill(pim_mom, pim_CTOF_vt);
							hi_pimc_vz.fill(pim_vz);
							hi_pim_CTOF_t.fill(pim_CTOF_t - STT - pim_vz/ (pimc_beta * 29.98f));
							hi_pim_CTOF_path.fill(pim_CTOF_path);
							
							
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
					if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut) {	//&& fdkpMass2Cut //&& research_fitter.dc_fiducial_cut(k, recbank, track_Bank, traj_Bank, run_Bank)
						if (recSCBank.getShort("pindex", r) == kp_part_ind && recSCBank.getByte("layer", r) == 2) {
							kp_FTOF_pad1b = recSCBank.getShort("component", r);
							kp_FTOF1b_t = recSCBank.getFloat("time", r);
							kp_FTOF1b_path = recSCBank.getFloat("path", r);
							float kp_beta = kp_mom / (float) Math.sqrt(kp_mom * kp_mom + PDGDatabase.getParticleById(321).mass() * PDGDatabase.getParticleById(321).mass());
							//kp_FTOF1b_vt = kp_FTOF1b_t - kp_FTOF1b_path / (k_ftb_beta * 29.98f) - STT - kp_vz/ (kp_ftb_beta * 29.98f);
							//kp_FTOF1b_vt = kp_FTOF1b_t - kp_FTOF1b_path / (kp_beta * 29.98f) - STT - kp_vz/ (kp_beta * 29.98f);
							//kp_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							float kp_FTOF1b_TOF = (float) kp_FTOF1b_t - STT - kp_vz/ (kp_beta * 29.98f);
							//difference in the measured and computed vertex times 
							kp_FTOF1b_vt = (float) kp_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass())));
							boolean fdkpvzCut = kp_vz > -10 && kp_vz < 2;
							float fkp_vtime_cut = Math.abs(0.05+Math.exp(-1.0*mom+0.4/mom));
							boolean fdkpCuts = (Math.abs(kp_FTOF1b_vt) <  fkp_vtime_cut && kp_FTOF1b_TOF > 20 && kp_FTOF1b_TOF < 55 && beta > 0.4 && beta < 1.05 && mom > 0.4 && mom < Eb && fdkpvzCut);
							if ( fdkpCuts ) {// Math.abs(kp_FTOF1b_vt) < 0.5 && kp_mom < 2.8 //fdkpvzCut && 
								Vkp = new Particle(pid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);

								if (runType == "mc"){
									//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
									//Random rand = new Random();
    								//double smearFactor = rand.nextGaussian()/100;
    								double smearedP = Vkp.p()*(1+smearFactor);
									Vkp.setP(smearedP);

									//Giovanni smearing function
									//LorentzVector lv_smeared = new LorentzVector();
									//lv_smeared = smear(Vkp.vector());
									//Vkp.setVector(lv_smeared, Vkp.vertex());
								}

								fkps.add(Vkp);
								kps.add(Vkp);
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
							boolean cdkpvzCut = kp_vz > -8 && kp_vz < 1;
							float kp_CTOF_TOF = (float) kp_CTOF_t - STT - kp_vz/ (kpc_beta * 29.98f);
							//difference in the measured and computed vertex times 
							kp_CTOF_vt = (float) kp_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(321).mass()*PDGDatabase.getParticleById(321).mass())));
							float ckp_vtime_ucut = 0.05+Math.exp(-6.0*mom+0.15/mom);
							float ckp_vtime_lcut = -(0.05+Math.exp(-3.5*mom+0.15/mom));
							boolean cdkpCut = kp_CTOF_vt < ckp_vtime_ucut && kp_CTOF_vt > ckp_vtime_lcut && kp_CTOF_TOF > 0.5 && kp_CTOF_TOF < 4 && beta > 0.2 && beta < 1.05 && mom > 0.2 && mom < 3.0 && cdkpvzCut;
							if (cdkpCut) { //cdkpvzCut && Math.abs(kp_CTOF_vt) < 0.4 && beta > 0.2 && beta < 1.05 && mom > 0.2 && mom < 3.0
								Vkpc = new Particle(pid, kp_px, kp_py, kp_pz, kp_vx, kp_vy, kp_vz);
								kps.add(Vkpc);
								ckps.add(Vkpc);
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
					if (recSCBank.getByte("detector", r) == 12 && inDC && fdChi2pidCut ) {	// && fdChi2pidCut //&& research_fitter.dc_fiducial_cut(k, recbank, track_Bank, traj_Bank, run_Bank)
						if (recSCBank.getShort("pindex", r) == km_part_ind && recSCBank.getByte("layer", r) == 2) {
							km_FTOF_pad1b = recSCBank.getShort("component", r);
							km_FTOF1b_t = recSCBank.getFloat("time", r);
							km_FTOF1b_path = recSCBank.getFloat("path", r);
							float km_beta = km_mom / (float) Math.sqrt(km_mom * km_mom + PDGDatabase.getParticleById(-321).mass() * PDGDatabase.getParticleById(-321).mass());
							//km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_ftb_beta * 29.98f) - STT - km_vz/ (km_ftb_beta * 29.98f);
							//km_FTOF1b_vt = km_FTOF1b_t - km_FTOF1b_path / (km_beta * 29.98f) - STT - km_vz/ (km_beta * 29.98f);
							//km_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							boolean fdkmvzCut = km_vz > -10 && km_vz < 2;
							float km_FTOF1b_TOF = (float) km_FTOF1b_t - STT - km_vz/ (km_beta * 29.98f);
							boolean km_FTOF1b_TOFCut = 22 < km_FTOF1b_TOF && km_FTOF1b_TOF < 28;
							//difference in the measured and computed vertex times 
							km_FTOF1b_vt = (float) km_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass())));
							float fkm_vtime_cut = Math.abs(0.05+Math.exp(-1.0*mom+0.4/mom));
							boolean fdkmCut = Math.abs(km_FTOF1b_vt) < fkm_vtime_cut && km_FTOF1b_TOF > 20 && km_FTOF1b_TOF < 35 && beta > 0.4 && beta < 1.05 && mom > 0.4 && mom < Eb && fdkmvzCut;							
							if (fdkmCut) {//&& fdkmvzCut && km_FTOF1b_TOFCut && Math.abs(km_FTOF1b_vt) < 0.5 && beta > 0.4 && beta < 1.05 && mom > 0.4 && mom < Eb
								Vkm = new Particle(pid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);
								
								if (runType == "mc"){
									//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
									//Random rand = new Random();
    								//double smearFactor = rand.nextGaussian()/100;
    								double smearedP = Vkm.p()*(1+smearFactor);
									Vkm.setP(smearedP);

									//Giovanni smearing function
									//LorentzVector lv_smeared = new LorentzVector();
									//lv_smeared = smear(Vkm.vector());
									//Vkm.setVector(lv_smeared, Vkm.vertex());

								}

								fkms.add(Vkm);
								kms.add(Vkm);
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
							boolean cdkmvzCut = kp_vz > -8 && kp_vz < 1;
							float km_CTOF_TOF = (float) km_CTOF_t - STT - km_vz/ (kmc_beta * 29.98f);
							//difference in the measured and computed vertex times 
							km_CTOF_vt = (float) km_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(-321).mass()*PDGDatabase.getParticleById(-321).mass())));
							float ckm_vtime_ucut = 0.05+Math.exp(-6.0*mom+0.15/mom);
							float ckm_vtime_lcut = -(0.05+Math.exp(-3.5*mom+0.15/mom));
							boolean cdkmCut = km_CTOF_vt < ckm_vtime_ucut && km_CTOF_vt > ckm_vtime_lcut && km_CTOF_TOF > 0.5 && km_CTOF_TOF < 4 && beta > 0.2 && beta < 1.05 && mom > 0.2 && mom < 3.0 && cdkmvzCut;
							if (cdkmCut) { //&& cdkmvzCut && Math.abs(km_CTOF_vt) < 0.4 && beta > 0.2 && beta < 1.05 && mom > 0.2 && mom < 3.0
								Vkmc = new Particle(pid, km_px, km_py, km_pz, km_vx, km_vy, km_vz);		
								kms.add(Vkmc);
								ckms.add(Vkmc);
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
							//prot_FTOF1b_vt = prot_FTOF1b_t - prot_FTOF1b_path / (prot_beta * 29.98f) - STT - prot_vz/ (prot_beta * 29.98f);
							//prot_FTOF1b_vt = STT - recFTbank.getFloat("vt", k);
							float prot_FTOF1b_TOF = (float) prot_FTOF1b_t - STT - prot_vz/ (prot_beta * 29.98f);
							boolean prot_FTOF1b_TOFCut = 22 < prot_FTOF1b_TOF && prot_FTOF1b_TOF < 32;
							//difference in the measured and computed vertex times 
							prot_FTOF1b_vt = (float) prot_FTOF1b_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass())));
							if (Math.abs(prot_FTOF1b_vt) < 0.5 && beta > 0.4 && beta < 1.05 && mom > 0.4 && mom < Eb) {// && prot_FTOF1b_TOFCut
								//Vprot = new LorentzVector(prot_px * proCor, prot_py* proCor, prot_pz* proCor, Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass()));
								Vprot = new Particle(pid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);
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
							float prot_CTOF_TOF = (float) prot_CTOF_t - STT - prot_vz/ (protc_beta * 29.98f);
							//difference in the measured and computed vertex times 
							prot_CTOF_vt = (float) prot_CTOF_TOF * (1 - Math.sqrt((mom*mom + mass)/(mom*mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass())));
							if (Math.abs(prot_CTOF_vt) < 0.4 && beta > 0.2 && beta < 1.05 && mom > 0.2 && mom < 3.0) {
								//Vprotc = new LorentzVector(prot_px * proCor, prot_py* proCor, prot_pz* proCor, Math.sqrt(prot_mom * prot_mom + PDGDatabase.getParticleById(2212).mass()*PDGDatabase.getParticleById(2212).mass()));
								Vprotc = new Particle(pid, prot_px, prot_py, prot_pz, prot_vx, prot_vy, prot_vz);		
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
			
			
			
			if ( q > 0 && ((pid == 45 && vz > -9 && vz < 6) || (pid == 2212 && vz > -9 && vz < 6) || (pid == 321 && vz > -9 && vz < 6) || (pid == 211 && vz > -9 && vz < 6) ) ) {//|| (pid == -11 && vz > -9 && vz < 6)
				if (inDC && fdChi2pidCut && mom > 0.6) {
					hi_fd_pos_mass.fill(mass);
					hi_FD_pos_beta_mom.fill(mom, beta);
					hi_FD_pos_mass_mom.fill(mom, mass);
					hi_FD_pos_mass_the.fill(the, mass);
					hi_FD_pos_mass_phi.fill(phi, mass);		
				}
				if (inCD && cdChi2pidCut && mom > 0.2) {
					hi_cd_pos_mass.fill(mass);
					hi_CD_pos_beta_mom.fill(mom, beta);
					hi_CD_pos_mass_mom.fill(mom, mass);
					hi_CD_pos_mass_the.fill(the, mass);
					hi_CD_pos_mass_phi.fill(phi, mass);

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
			if ( q < 0 && ((pid == -2212 && vz > -9 && vz < 6) || (pid == -321 && vz > -9 && vz < 6) || (pid == -211 && vz > -9 && vz < 6))) {//|| || (pid == -11 && vz > -9 && vz < 6)
				if (inDC && fdChi2pidCut && mom > 0.6) {
					hi_fd_neg_mass.fill(mass);
					hi_FD_neg_beta_mom.fill(mom, beta);
					hi_FD_neg_mass_mom.fill(mom, mass);
					hi_FD_neg_mass_the.fill(the, mass);
					hi_FD_neg_mass_phi.fill(phi, mass);		
				}
				if (inCD && cdChi2pidCut && mom > 0.2) {
					hi_cd_neg_mass.fill(mass);
					hi_CD_neg_beta_mom.fill(mom, beta);
					hi_CD_neg_mass_mom.fill(mom, mass);
					hi_CD_neg_mass_the.fill(the, mass);
					hi_CD_neg_mass_phi.fill(phi, mass);				
				}

			}
			if (q == 0 ) {
				if (inDC && fdChi2pidCut) {
					hi_FD_neutral_beta_mom.fill(mom, beta);
					hi_FD_neutral_mass_mom.fill(mom, mass);
					hi_FD_neutral_mass_the.fill(the, mass);
					hi_FD_neg_mass_phi.fill(phi, mass);		
				}
				if (inCD && cdChi2pidCut) {
					hi_CD_neutral_beta_mom.fill(mom, beta);
					hi_CD_neutral_mass_mom.fill(mom, mass);
					hi_CD_neutral_mass_the.fill(the, mass);
					hi_CD_neutral_mass_phi.fill(phi, mass);				
				}
			}
			
		} // FOR LOOP

	
	} //MAKEOTHER
	

	public boolean select_efkpfkpfkm(){
		boolean res = false;
		efkpfkpfkm_found_lambda = false;
		if(found_eFD && fkps.size() >= 2 && fkms.size() >= 1 ){
			//double avg_vz = (double)(fkps.get(0).vz() + fkps.get(1).vz() + fkms.get(0).vz())/3;
			//Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(fkps.get(1).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			efkpfkpfkm_MM_efkpfkp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Ve);Vmissekpkpkm.sub(fkps.get(0).vector());Vmissekpkpkm.sub(fkps.get(1).vector());Vmissekpkpkm.sub(fkms.get(0).vector());
			efkpfkpfkm_MM_efkpfkpfkm = (float)Vmissekpkpkm.mass();

			if (Math.abs(efkpfkpfkm_MM_efkpfkpfkm - 1.136) <= 0.093){ //Math.abs(efkpfkpfkm_MM_efkpfkpfkm - 1.132) <= 0.0766 (3sigma cut, 3*0.0255)//efkpfkpfkm_MM_efkpfkpfkm < 1.17 && efkpfkpfkm_MM_efkpfkpfkm > 1.06(to selectLambdaK-decay only/excluds 2030)//efkpfkpfkm_MM_efkpfkpfkm < 1.28 && efkpfkpfkm_MM_efkpfkpfkm > 1.06(to include SigmaK- decay from 2030 as well)// efkpfkpfkm_MM_efkpfkpfkm < 1.27 && efkpfkpfkm_MM_efkpfkpfkm > 1.05 //Math.abs(efkpfkpfkm_MM_efkpfkpfkm - 1.115683) <= 0.052
				efkpfkpfkm_found_lambda = true;
				//constraint Lambda with nominal mass
				LorentzVector Vconstrainedlambda = new LorentzVector(Vmissekpkpkm.px(), Vmissekpkpkm.py(), Vmissekpkpkm.pz(), Math.sqrt(Vmissekpkpkm.p()*Vmissekpkpkm.p() + 1.115683f*1.115683f));
				LorentzVector VxifromLambda = new LorentzVector(0, 0, 0, 0);
				VxifromLambda.add(Vconstrainedlambda);VxifromLambda.add(fkms.get(0).vector());
				//VxifromLambda.add(Vmissekpkpkm);VxifromLambda.add(fkms.get(0).vector());
				
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
		if(found_eFD && fkps.size() >= 1&& ckps.size()>=1 && fkms.size() >= 1 ){

			LorentzVector lv_missekpkpkm_nocorn = new LorentzVector(0, 0, 0, 0);
			lv_missekpkpkm_nocorn.add(VT);lv_missekpkpkm_nocorn.add(VB);lv_missekpkpkm_nocorn.sub(Ve);lv_missekpkpkm_nocorn.sub(fkps.get(0).vector());lv_missekpkpkm_nocorn.sub(ckps.get(0).vector());lv_missekpkpkm_nocorn.sub(fkms.get(0).vector());
			efkpckpfkm_MM_efkpckpfkm_nocorn = (float)lv_missekpkpkm_nocorn.mass();
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

			//double avg_vz = (double)(fkps.get(0).vz() + ckps.get(0).vz() + fkms.get(0).vz())/3;
			//Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);

			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			efkpckpfkm_MM_efkpckp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Ve);Vmissekpkpkm.sub(fkps.get(0).vector());Vmissekpkpkm.sub(ckps.get(0).vector());Vmissekpkpkm.sub(fkms.get(0).vector());
			efkpckpfkm_MM_efkpckpfkm = (float)Vmissekpkpkm.mass();

			if (Math.abs(efkpckpfkm_MM_efkpckpfkm - 1.115683) <= 0.052){
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
		if(found_eFD && ckps.size()>=2 && fkms.size() >= 1 ){
			//double avg_vz = (double)(ckps.get(0).vz() + ckps.get(1).vz() + fkms.get(0).vz())/3;
			//Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(ckps.get(0).vector());Vmissekpkp.sub(ckps.get(1).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			eckpckpfkm_MM_eckpckp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Ve);Vmissekpkpkm.sub(ckps.get(0).vector());Vmissekpkpkm.sub(ckps.get(1).vector());Vmissekpkpkm.sub(fkms.get(0).vector());
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
		if(found_eFD && fkps.size() >= 2 && ckms.size() >= 1 ){
			//double avg_vz = (double)(fkps.get(0).vz() + fkps.get(1).vz() + ckms.get(0).vz())/3;
			//Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(fkps.get(1).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			efkpfkpckm_MM_efkpfkp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Ve);Vmissekpkpkm.sub(fkps.get(0).vector());Vmissekpkpkm.sub(fkps.get(1).vector());Vmissekpkpkm.sub(ckms.get(0).vector());
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
		if(found_eFD && fkps.size() >= 1&& ckps.size()>=1 && ckms.size() >= 1 ){
			//double avg_vz = (double)(fkps.get(0).vz() + ckps.get(0).vz() + ckms.get(0).vz())/3;
			//Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(fkps.get(0).vector());Vmissekpkp.sub(ckps.get(0).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			efkpckpckm_MM_efkpckp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Ve);Vmissekpkpkm.sub(fkps.get(0).vector());Vmissekpkpkm.sub(ckps.get(0).vector());Vmissekpkpkm.sub(ckms.get(0).vector());
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
		if(found_eFD && ckps.size() >= 2 && ckms.size() >= 1 ){
			//double avg_vz = (double)(ckps.get(0).vz() + ckps.get(1).vz() + ckms.get(0).vz())/3;
			//Vftel_corrected = electron_corn(ftels.get(0), e_ftCal_hitPosition, avg_vz);
			// with the common vertex correction applied to electron
			LorentzVector Vmissekpkp = new LorentzVector(0, 0, 0, 0);
			Vmissekpkp.add(VT);Vmissekpkp.add(VB);Vmissekpkp.sub(Ve);Vmissekpkp.sub(ckps.get(0).vector());Vmissekpkp.sub(ckps.get(1).vector());
			//ekpkp_MM_req_km = (float)Vmissekpkp.mass();
			eckpckpckm_MM_eckpckp = (float)Vmissekpkp.mass();

			LorentzVector Vmissekpkpkm = new LorentzVector(0, 0, 0, 0);
			Vmissekpkpkm.add(VT);Vmissekpkpkm.add(VB);Vmissekpkpkm.sub(Ve);Vmissekpkpkm.sub(ckps.get(0).vector());Vmissekpkpkm.sub(ckps.get(1).vector());Vmissekpkpkm.sub(ckms.get(0).vector());
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
		found_fastpip = false;
		found_slowpip = false;
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
		*/
			
	}

	public void processEvent(DataEvent event) {
		resetCounters();
		//Timothy's analysis fitter
        //GenericKinematicFitter fitter = new analysis_fitter(Eb);

  //if (research_fitter.banks_test(event)) {
		//if (event.hasBank("RECFT::Event"))
		//	fillRecBank(event.getBank("RECFT::Event"));
		if (event.hasBank("REC::Event"))
			fillRecBank(event.getBank("REC::Event"));
		//if (event.hasBank("REC::Particle") && event.hasBank("RECFT::Particle") && event.hasBank("REC::ForwardTagger")) e_ft_part_ind = makeFTElectron(event.getBank("REC::Particle"), event.getBank("RECFT::Particle"), event.getBank("REC::ForwardTagger"));
		
		//if (event.hasBank("REC::Particle")) e_ft_part_ind = makeFDElectron(event.getBank("REC::Particle"));

		// load the hipo banks
        // assumption is we are using trains which would require all of these banks to exist
        HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle"); // load particle bank
        HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
        HipoDataBank cc_Bank = (HipoDataBank) event.getBank("REC::Cherenkov");
        HipoDataBank track_Bank = (HipoDataBank) event.getBank("REC::Track");
        HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");
        HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config");
        HipoDataBank scint_Bank = (HipoDataBank) event.getBank("REC::Scintillator");

		if (event.hasBank("REC::Particle")) e_ft_part_ind = makeFDElectron(rec_Bank, cal_Bank, track_Bank, traj_Bank, run_Bank, cc_Bank);
		if (event.hasBank("MC::Particle") == true) fillMCPartBank(event.getBank("MC::Particle"), event.getBank("RUN::config"));
		if(e_ft_part_ind > -1 && found_eFD ) {	//&& scint_Bank.getByte("sector", e_ft_part_ind) == 2	
			//if (event.hasBank("REC::Particle") && event.hasBank("REC::Scintillator")) makeOthers(event.getBank("REC::Particle"), event.getBank("REC::Scintillator"));
			//HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank, HipoDataBank cc_Bank
			if (event.hasBank("REC::Particle") && event.hasBank("REC::Scintillator")) makeOthers(event.getBank("REC::Particle"), event.getBank("REC::Scintillator"), cal_Bank, track_Bank, traj_Bank, run_Bank, cc_Bank);
			//if (event.hasBank("REC::Scintillator")) fillFTOF(event.getBank("REC::Particle"), event.getBank("REC::Scintillator"));
			//System.out.println("found pip :: " + pips.size() + " pips tracks in FD & " + cpips.size() + " pip tracks in CD.");
			FillHists();
			
		} //e_ft_part_ind > -1
  //}
		
	} //processEvent


	/*
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
    //	fitxi(dg_mm2.getH1F("hi_"))


	}
	*/


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

			/*
			H_FT_e_t_f.fill(Math.toDegrees(Ve.phi()), Math.toDegrees(Ve.theta()));
			H_FT_e_p_f.fill(Math.toDegrees(Ve.phi()), Ve.p());
			H_FT_e_p_the.fill(Math.toDegrees(Ve.theta()), Ve.p());
			hi_rec_fde_vz.fill(fdels.get(0).vz());
			hi_rec_fde_theta_vz.fill(fdels.get(0).vz(), Math.toDegrees(fdels.get(0).theta()));
			

//			H_FT_e_t_f.fill(e_phi, e_the);
//			H_FT_e_p_f.fill(e_phi, e_mom);
//			H_FT_e_p_the.fill(e_the, e_mom);
			H_FT_W_Q2.fill(e_W, e_Q2);
			H_FT_W.fill(e_W);
			H_FT_Q2.fill(e_Q2);	
			H_virphoton.fill(e_virphoton);	
			//*/


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
	

		// ekpkp[km required 
		if(select_efkpfkpfkm() && !select_efkpckpfkm() && !select_efkpfkpckm()){
			
			
			hi_efkpfkpfkm_MM_efkpfkp.fill(efkpfkpfkm_MM_efkpfkp);
			hi_efkpfkpfkm_MM_efkpfkpfkm.fill(efkpfkpfkm_MM_efkpfkpfkm);
			hi_efkpfkpfkm_MM_efkpfkp_MM_efkpfkpfkm.fill(efkpfkpfkm_MM_efkpfkpfkm, efkpfkpfkm_MM_efkpfkp);

			if(efkpfkpfkm_found_lambda) {
				hi_efkpfkpfkm_MM_efkpfkp_lam_evnt.fill(efkpfkpfkm_MM_efkpfkp);
				hi_efkpfkpfkm_IM_kmlambda.fill(efkpfkpfkm_IM_kmlambda);
			}


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


			// set up boost to gamma*-nucleon center of mass frame
			LorentzVector lv_q = new LorentzVector(VB); lv_q.sub(Ve);
        	LorentzVector gN = new LorentzVector(lv_q);
			gN.add(VT);
			Vector3 gNBoost = gN.boostVector();
			gNBoost.negative();

			// boost to gamma*-nucleon center of mass frame
			LorentzVector lv_target_gN = new LorentzVector(VT); lv_target_gN.boost(gNBoost);
			LorentzVector lv_e_gN = new LorentzVector(Ve); lv_e_gN.boost(gNBoost);
			LorentzVector lv_fastkp_gN = new LorentzVector(Vfastkp.vector()); lv_fastkp_gN.boost(gNBoost);
			LorentzVector lv_slowkp_gN = new LorentzVector(Vslowkp.vector()); lv_slowkp_gN.boost(gNBoost);
			LorentzVector lv_km_gN = new LorentzVector(fkms.get(0).vector()); lv_km_gN.boost(gNBoost);


			// two kps dtheta check in CM frame
			///*
			hi_rec_kps_dtheta.fill(Math.toDegrees(lv_slowkp_gN.theta()) - Math.toDegrees(lv_fastkp_gN.theta()));
			hi_rec_kps_dtheta_slowkp_theta.fill(Math.toDegrees(lv_slowkp_gN.theta()), Math.toDegrees(lv_slowkp_gN.theta()) - Math.toDegrees(lv_fastkp_gN.theta()));
			//hi_rec_kps_dtheta_fastkp_theta.fill(Math.toDegrees(Vfastkp.theta()), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			//hi_rec_kps_dtheta_slowkp_phi.fill(Math.toDegrees(Vslowkp.phi()), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			//hi_rec_kps_dtheta_fastkp_phi.fill(Math.toDegrees(Vfastkp.phi()), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			//hi_rec_kps_dtheta_slowkp_p.fill(Vslowkp.p(), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			//hi_rec_kps_dtheta_fastkp_p.fill(Vfastkp.p(), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			//*/

			// two kps dtheta check in Lab frame
			/*
			//hi_rec_kps_dtheta.fill(Math.abs(Math.toDegrees(Vfastkp.theta()) - Math.toDegrees(Vslowkp.theta())));
			hi_rec_kps_dtheta_slowkp_theta.fill(Math.toDegrees(Vslowkp.theta()), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			hi_rec_kps_dtheta_fastkp_theta.fill(Math.toDegrees(Vfastkp.theta()), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			hi_rec_kps_dtheta_slowkp_phi.fill(Math.toDegrees(Vslowkp.phi()), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			hi_rec_kps_dtheta_fastkp_phi.fill(Math.toDegrees(Vfastkp.phi()), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			hi_rec_kps_dtheta_slowkp_p.fill(Vslowkp.p(), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			hi_rec_kps_dtheta_fastkp_p.fill(Vfastkp.p(), Math.toDegrees(Vslowkp.theta()) - Math.toDegrees(Vfastkp.theta()));
			//*/

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


			/*
			//Store three kaon event information in the arraylist so that we can implement event mixing technique later on
            //public static ArrayList<LorentzVector> lv_my_fdels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
            lv_my_fdels.add(Ve);
            lv_my_fast_fkps.add(Vfastkp.vector());
            lv_my_slow_fkps.add(Vslowkp.vector());
            lv_my_fkms.add(fkms.get(0).vector());
			//*/

			///*
			//public static ArrayList<LorentzVector> lv_my_xi_fdels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms;
			//Xi data sample from MM(eK+K+K-) cut to select lambda/sigma region from the fit and store information so tha we can implement event mixing technique later on to estimate bg shape
			//1.11 <= mm_fkpskpkm && mm_fkpskpkm <= 1.21
			/*
			if (1.05 <= mm_fkpskpkm && mm_fkpskpkm <= 1.29){ //f18out:: 1.07 <= mm_fkpskpkm && mm_fkpskpkm <= 1.25 // f18in:: 1.02 <= mm_fkpskpkm && mm_fkpskpkm <= 1.28// s19:: 1.05 <= mm_fkpskpkm && mm_fkpskpkm <= 1.29
				lv_my_xi_fdels.add(Ve);
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
			hi_efkpfkpfkm_MM_efastkpaspipslowkpfkm.fill(lv_mm_fkpaspipskpkm.mass(), efkpfkpfkm_MM_efkpfkpfkm);
			double mm_fkpaspipskpkm = lv_mm_fkpaspipskpkm.mass();

			//// considering slow kp as pip
			LorentzVector lv_mm_fkpskpaspipkm = new LorentzVector();lv_mm_fkpskpaspipkm.add(lv_mm_fkpkm);lv_mm_fkpskpaspipkm.sub(lv_skpaspip);
			hi_efkpfkpfkm_MM_eslowkpaspipfastkpfkm.fill(lv_mm_fkpskpaspipkm.mass(), efkpfkpfkm_MM_efkpfkpfkm);
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
			double nu = VB.e()-Ve.e();
			double x  = q2 / (2 * PDGDatabase.getParticleById(2212).mass() * nu);
			double W  = Math.pow(Math.pow(PDGDatabase.getParticleById(2212).mass(),2)+2*PDGDatabase.getParticleById(2212).mass()*nu - q2, 0.5);

			// lab e kinematics
			double beam_e = VB.e();
			double e_px = Ve.px();
			double e_py = Ve.py();
			double e_pz = Ve.pz();
			double e_p = Ve.p();
			double e_e = Ve.e();
			double e_vx = Vfde.vx();
			double e_vy = Vfde.vy();
			double e_vz = Vfde.vz();
			double e_theta = Math.toDegrees(Ve.theta());
			double e_phi = Math.toDegrees(Ve.phi());
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


			/* // to save all three kaon events

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
			//*/


		//	/* // to save Lambda/Sigma0 events in the missing mass

			//Calculate Lambda/Sigma/BG probability for each event in the selected region (using MM(KKK)) according to fitting functions

			double global = fn_lam_fit.evaluate(mm_fkpskpkm) + fn_sig_fit.evaluate(mm_fkpskpkm) + fn_bg_fit.evaluate(mm_fkpskpkm);
			double w_lamsig = fn_lamsig_fit.evaluate(mm_fkpskpkm)/global;
			double w_lam = fn_lam_fit.evaluate(mm_fkpskpkm)/global;
			double w_sig = fn_sig_fit.evaluate(mm_fkpskpkm)/global;
			double w_bg = fn_bg_fit.evaluate(mm_fkpskpkm)/global;

			outFile.append(runnum+" "+evnum+" "+beam_e+" "+q2+" "+nu+" "+x+" "+W+" ");
			outFile.append(e_px+" "+e_py+" "+e_pz+" "+e_p+" "+e_e+" "+e_vx+" "+e_vy+" "+e_vz+" "+e_theta+" "+e_phi+" ");
			outFile.append(fkp_px+" "+fkp_py+" "+fkp_pz+" "+fkp_p+" "+fkp_e+" "+fkp_vx+" "+fkp_vy+" "+fkp_vz+" "+fkp_theta+" "+fkp_phi+" ");
			outFile.append(skp_px+" "+skp_py+" "+skp_pz+" "+skp_p+" "+skp_e+" "+skp_vx+" "+skp_vy+" "+skp_vz+" "+skp_theta+" "+skp_phi+" ");
			outFile.append(km_px+" "+km_py+" "+km_pz+" "+km_p+" "+km_e+" "+km_vx+" "+km_vy+" "+km_vz+" "+km_theta+" "+km_phi+" ");
			outFile.append(m_fkpskp+" "+m_fkpkm+" "+m_skpkm+" "+m_fkpskpkm+" ");
			outFile.append(mm_e+" "+mm_fkp+" "+mm_skp+" "+mm_km+" "+mm_fkpskp+" "+mm_fkpkm+" "+mm_skpkm+" "+mm_fkpskpkm+" ");
			outFile.append(mm_fkpaspipskpkm+" "+mm_fkpskpaspipkm+" "+mm_fkpskpkmaspim+" "+mm_fkpaspipskpaspipkm+" "+mm_fkpaspipskpkmaspim+" "+mm_fkpskpaspipkmaspim+" "+mm_fkpaspipskpaspipkmaspim+" "+w_lam+" "+w_sig+" "+w_bg+" "+w_lamsig+"\n");	
			//System.out.println(w_lam + " " + w_sig + " " + w_bg);
			//*/
			//0.16 <=q2 && q2 <= 1.28 (f18outfde), 1.28 <=q2 && q2 <= 2.88 (other)

			if (1.05<=mm_fkpskpkm && mm_fkpskpkm<=1.29 && Math.toDegrees(Ve.theta()) >= 5.0 && Math.toDegrees(Ve.theta()) <= 35.0 && Ve.p() >= 1.0 && Math.toDegrees(Vfastkp.theta()) <= 35.0 && Math.toDegrees(Vslowkp.theta()) <= 35.0 && Math.toDegrees(fkms.get(0).theta()) <= 35.0 && Math.toDegrees(Ve.theta()) <= 35.0  && Math.toDegrees(Vfastkp.theta()) >= 5.0 && Math.toDegrees(Vslowkp.theta()) >= 5.0 && Math.toDegrees(fkms.get(0).theta()) >= 5.0){
				//*/
				H_FT_e_t_f.fill(Math.toDegrees(Ve.phi()), Math.toDegrees(Ve.theta()));
				H_FT_e_p_f.fill(Math.toDegrees(Ve.phi()), Ve.p());
				H_FT_e_p_the.fill(Math.toDegrees(Ve.theta()), Ve.p());
				hi_rec_fde_vz.fill(fdels.get(0).vz());
				hi_rec_fde_theta_vz.fill(fdels.get(0).vz(), Math.toDegrees(fdels.get(0).theta()));
			
				H_FT_W_Q2.fill(e_W, e_Q2);
				H_FT_W.fill(e_W);
				H_FT_Q2.fill(e_Q2);	
				H_virphoton.fill(e_virphoton);
				//*/
			
				hi_rec_kp1_p_the.fill(Math.toDegrees(Vfastkp.theta()), Vfastkp.p());
				hi_rec_kp1_p_phi.fill(Math.toDegrees(Vfastkp.phi()), Vfastkp.p());
				hi_rec_kp1_the_phi.fill(Math.toDegrees(Vfastkp.phi()), Math.toDegrees(Vfastkp.theta()));
				hi_rec_kp2_p_the.fill(Math.toDegrees(Vslowkp.theta()), Vslowkp.p());
				hi_rec_kp2_p_phi.fill(Math.toDegrees(Vslowkp.phi()), Vslowkp.p());
				hi_rec_kp2_the_phi.fill(Math.toDegrees(Vslowkp.phi()), Math.toDegrees(Vslowkp.theta()));
				hi_rec_km_p_the.fill(Math.toDegrees(fkms.get(0).theta()), fkms.get(0).p());
				hi_rec_km_p_phi.fill(Math.toDegrees(fkms.get(0).phi()), fkms.get(0).p());
				hi_rec_km_the_phi.fill(Math.toDegrees(fkms.get(0).phi()), Math.toDegrees(fkms.get(0).theta()));
				hi_rec_kp1_vz.fill(Vfastkp.vz());
				hi_rec_kp2_vz.fill(Vslowkp.vz());
				hi_rec_km_vz.fill(fkms.get(0).vz());
				hi_rec_kp1_theta_vz.fill(Vfastkp.vz(), Math.toDegrees(Vfastkp.theta()));
				hi_rec_kp2_theta_vz.fill(Vslowkp.vz(), Math.toDegrees(Vslowkp.theta()));
				hi_rec_km_theta_vz.fill(fkms.get(0).vz(), Math.toDegrees(fkms.get(0).theta()));

			}

			if (Math.toDegrees(Ve.theta()) >= 5.0 && Math.toDegrees(Ve.theta()) <= 35.0 && Ve.p() >= 1.0 && Math.toDegrees(Vfastkp.theta()) <= 35.0 && Math.toDegrees(Vslowkp.theta()) <= 35.0 && Math.toDegrees(fkms.get(0).theta()) <= 35.0 && Math.toDegrees(Ve.theta()) <= 35.0  && Math.toDegrees(Vfastkp.theta()) >= 5.0 && Math.toDegrees(Vslowkp.theta()) >= 5.0 && Math.toDegrees(fkms.get(0).theta()) >= 5.0) {
				//	/*
				//Store three kaon event information in the arraylist so that we can implement event mixing technique later on
            	//public static ArrayList<LorentzVector> lv_my_fdels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
            	lv_my_fdels.add(Ve);
            	lv_my_fast_fkps.add(Vfastkp.vector());
            	lv_my_slow_fkps.add(Vslowkp.vector());
            	lv_my_fkms.add(fkms.get(0).vector());
				//*/
			}

			


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

				// labeling reconstructed kps with momentum
				Particle Vslowkp;
				Particle Vfastkp;

				if (fkps.get(0).p() > ckps.get(0).p()){
					Vfastkp = new Particle(fkps.get(0));
					Vslowkp = new Particle(ckps.get(0));
				} else {
					Vfastkp = new Particle(ckps.get(0));
					Vslowkp = new Particle(fkps.get(0));
				}

				/*
				hi_rec_kp1_p_the.fill(Math.toDegrees(Vfastkp.theta()), Vfastkp.p());
				hi_rec_kp1_p_phi.fill(Math.toDegrees(Vfastkp.phi()), Vfastkp.p());
				hi_rec_kp1_the_phi.fill(Math.toDegrees(Vfastkp.phi()), Math.toDegrees(Vfastkp.theta()));
				hi_rec_kp2_p_the.fill(Math.toDegrees(Vslowkp.theta()), Vslowkp.p());
				hi_rec_kp2_p_phi.fill(Math.toDegrees(Vslowkp.phi()), Vslowkp.p());
				hi_rec_kp2_the_phi.fill(Math.toDegrees(Vslowkp.phi()), Math.toDegrees(Vslowkp.theta()));
				hi_rec_km_p_the.fill(Math.toDegrees(fkms.get(0).theta()), fkms.get(0).p());
				hi_rec_km_p_phi.fill(Math.toDegrees(fkms.get(0).phi()), fkms.get(0).p());
				hi_rec_km_the_phi.fill(Math.toDegrees(fkms.get(0).phi()), Math.toDegrees(fkms.get(0).theta()));
				hi_rec_kp1_vz.fill(Vfastkp.vz());
				hi_rec_kp2_vz.fill(Vslowkp.vz());
				hi_rec_km_vz.fill(fkms.get(0).vz());
				hi_rec_kp1_theta_vz.fill(Vfastkp.vz(), Math.toDegrees(Vfastkp.theta()));
				hi_rec_kp2_theta_vz.fill(Vslowkp.vz(), Math.toDegrees(Vslowkp.theta()));
				hi_rec_km_theta_vz.fill(fkms.get(0).vz(), Math.toDegrees(fkms.get(0).theta()));
				*/
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
		if(select_efkpfkpckm() && !select_efkpckpfkm() && !select_efkpfkpfkm()){
			hi_efkpfkpckm_MM_efkpfkp.fill(efkpfkpckm_MM_efkpfkp);
			hi_efkpfkpckm_MM_efkpfkpckm.fill(efkpfkpckm_MM_efkpfkpckm);
			hi_efkpfkpckm_MM_efkpfkp_MM_efkpfkpckm.fill(efkpfkpckm_MM_efkpfkpckm, efkpfkpckm_MM_efkpfkp);
			if(efkpfkpckm_found_lambda){
				hi_efkpfkpckm_MM_efkpfkp_lam_evnt.fill(efkpfkpckm_MM_efkpfkp);
				hi_efkpfkpckm_IM_kmlambda.fill(efkpfkpckm_IM_kmlambda);
			}
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
			/*
			hi_rec_kp1_p_the.fill(Math.toDegrees(Vfastkp.theta()), Vfastkp.p());
			hi_rec_kp1_p_phi.fill(Math.toDegrees(Vfastkp.phi()), Vfastkp.p());
			hi_rec_kp1_the_phi.fill(Math.toDegrees(Vfastkp.phi()), Math.toDegrees(Vfastkp.theta()));
			hi_rec_kp2_p_the.fill(Math.toDegrees(Vslowkp.theta()), Vslowkp.p());
			hi_rec_kp2_p_phi.fill(Math.toDegrees(Vslowkp.phi()), Vslowkp.p());
			hi_rec_kp2_the_phi.fill(Math.toDegrees(Vslowkp.phi()), Math.toDegrees(Vslowkp.theta()));
			hi_rec_km_p_the.fill(Math.toDegrees(ckms.get(0).theta()), ckms.get(0).p());
			hi_rec_km_p_phi.fill(Math.toDegrees(ckms.get(0).phi()), ckms.get(0).p());
			hi_rec_km_the_phi.fill(Math.toDegrees(ckms.get(0).phi()), Math.toDegrees(ckms.get(0).theta()));
			hi_rec_kp1_vz.fill(Vfastkp.vz());
			hi_rec_kp2_vz.fill(Vslowkp.vz());
			hi_rec_km_vz.fill(ckms.get(0).vz());
			hi_rec_kp1_theta_vz.fill(Vfastkp.vz(), Math.toDegrees(Vfastkp.theta()));
			hi_rec_kp2_theta_vz.fill(Vslowkp.vz(), Math.toDegrees(Vslowkp.theta()));
			hi_rec_km_theta_vz.fill(ckms.get(0).vz(), Math.toDegrees(ckms.get(0).theta()));
			*/

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

		if(select_eckpckpckm()){
			hi_eckpckpckm_MM_eckpckp.fill(eckpckpckm_MM_eckpckp);
			hi_eckpckpckm_MM_eckpckpckm.fill(eckpckpckm_MM_eckpckpckm);
			hi_eckpckpckm_MM_eckpckp_MM_eckpckpckm.fill(eckpckpckm_MM_eckpckpckm, eckpckpckm_MM_eckpckp);
			if(eckpckpckm_found_lambda){
					hi_eckpckpckm_MM_eckpckp_lam_evnt.fill(eckpckpckm_MM_eckpckp);
					hi_eckpckpckm_IM_kmlambda.fill(eckpckpckm_IM_kmlambda);
			}
		}

			
	}
		
	public void plot(){

		myCanvas = new EmbeddedCanvasTabbed("Kinematics","REC-electron","REC-kp1(fast)","REC-kp2(slow)","REC-kps","REC-km","MM-ekpkpkm","MM-ekpkpkm1","MM-kpAspip","MM-scatter","MM-lam-evnt","IM-kmlambda","MM-spectra","FD-TOF", "CD-TOF", "CD-Part", "VTime", "particle-Vz","TOF-t","TOF-path", "Counter");
		

		myCanvas.getCanvas("Kinematics").divide(2, 2);
		myCanvas.getCanvas("Kinematics").setSize(1600, 1000);
		myCanvas.getCanvas("Kinematics").setGridX(false);
		myCanvas.getCanvas("Kinematics").setGridY(false);
		myCanvas.getCanvas("Kinematics").setAxisFontSize(18);
		myCanvas.getCanvas("Kinematics").setAxisTitleSize(24);
		myCanvas.getCanvas("Kinematics").draw(dg_kinematics);
		myCanvas.getCanvas("Kinematics").getPad(2).getAxisZ().setLog(true);
		
		//reconstructed electron

		myCanvas.getCanvas("REC-electron").divide(3, 2);
		myCanvas.getCanvas("REC-electron").setSize(1600, 1000);
		myCanvas.getCanvas("REC-electron").setGridX(false);
		myCanvas.getCanvas("REC-electron").setGridY(false);
		myCanvas.getCanvas("REC-electron").setAxisFontSize(18);
		myCanvas.getCanvas("REC-electron").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-electron").draw(dg_rec_electron);
		myCanvas.getCanvas("REC-electron").getPad(0).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(1).getAxisZ().setLog(true);
		myCanvas.getCanvas("REC-electron").getPad(2).getAxisZ().setLog(true);

		// reconstructed kp1
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
		myCanvas.getCanvas("REC-kp1(fast)").getPad(4).getAxisZ().setLog(true);
		//myCanvas.getCanvas("REC-kp1(fast)").getPad(10).getAxisZ().setLog(true);
		//myCanvas.getCanvas("REC-kp1(fast)").getPad(11).getAxisZ().setLog(true);

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
		myCanvas.getCanvas("REC-kp2(slow)").getPad(4).getAxisZ().setLog(true);
		//myCanvas.getCanvas("REC-kp2(slow)").getPad(10).getAxisZ().setLog(true);
		//myCanvas.getCanvas("REC-kp2(slow)").getPad(11).getAxisZ().setLog(true);


		// two K+s kps
		myCanvas.getCanvas("REC-kps").divide(3,4);
		myCanvas.getCanvas("REC-kps").setSize(1600, 1000);
		myCanvas.getCanvas("REC-kps").setGridX(false);
		myCanvas.getCanvas("REC-kps").setGridY(false);
		myCanvas.getCanvas("REC-kps").setAxisFontSize(18);
		myCanvas.getCanvas("REC-kps").setAxisTitleSize(24);
		myCanvas.getCanvas("REC-kps").draw(dg_rec_kps);
		

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
		myCanvas.getCanvas("REC-km").getPad(4).getAxisZ().setLog(true);
		//myCanvas.getCanvas("REC-km").getPad(10).getAxisZ().setLog(true);
		//myCanvas.getCanvas("REC-km").getPad(11).getAxisZ().setLog(true);

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

		myCanvas.getCanvas("MM-kpAspip").divide(3, 2);
		myCanvas.getCanvas("MM-kpAspip").setSize(1600, 1000);
		myCanvas.getCanvas("MM-kpAspip").setGridX(false);
		myCanvas.getCanvas("MM-kpAspip").setGridY(false);
		myCanvas.getCanvas("MM-kpAspip").setAxisFontSize(18);
		myCanvas.getCanvas("MM-kpAspip").setAxisTitleSize(24);
		myCanvas.getCanvas("MM-kpAspip").draw(dg_rec_kpaspip);

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


		// TOF
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

		//### CD Particle
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
				
	}

	
	public void save() {
		
		//TDirectory dirout = new TDirectory();
		//dirout.mkdir("");
		//dirout.cd("/FTElec/");
		//dirout.addDataSet(H_FT_e_t_f, H_FT_e_p_f, H_FT_W_Q2);
		myCanvas.getCanvas("Kinematics").save("kinematics.png");
		myCanvas.getCanvas("REC-electron").save("electron.png");
		myCanvas.getCanvas("REC-kp1(fast)").save("fast_kaon_plus.png");
		myCanvas.getCanvas("REC-kp2(slow)").save("slow_kaon_plus.png");
		myCanvas.getCanvas("REC-km").save("kaon_minus.png");
		//myCanvas.getCanvas("electron-res").save("ft_electron_resolution.png");
		//myCanvas.getCanvas("kp1(fast-kp)-res").save("fast_kaon_plus_resolution.png");
		//myCanvas.getCanvas("kp2(slow-kp)-res").save("slow_kaon_plus_resolution.png");
		//myCanvas.getCanvas("REC-Xi").save("cascade.png");
		myCanvas.getCanvas("MM-ekpkpkm").save("ekpkpkm_mm_ekpkp_phasespace.png");
		myCanvas.getCanvas("MM-ekpkpkm1").save("ekpkpkm_mm_ekpkpkm_phasespace.png");
		myCanvas.getCanvas("MM-scatter").save("ekpkpkm_mm_ekpkp_mm_ekpkpkm_scatter.png");
		myCanvas.getCanvas("MM-lam-evnt").save("ekpkpkm_mm_ekpkp_lambdaevent.png");
		myCanvas.getCanvas("IM-kmlambda").save("ekpkpkm_im_kmlambda.png")
		myCanvas.getCanvas("MM-spectra").save("mm_ekpkpkm.png");
		myCanvas.getCanvas("CD-TOF").save("cd_beta_p_overview.png");
		myCanvas.getCanvas("FD-TOF").save("fd_beta_p_overview.png");
		myCanvas.getCanvas("VTime").save("particle_dtime_overview.png");
		myCanvas.getCanvas("particle-Vz").save("particle_vz_overview.png");
		myCanvas.getCanvas("TOF-t").save("tof_time_overview.png");
		myCanvas.getCanvas("TOF-path").save("tof_path_overview.png");
		myCanvas.getCanvas("Counter").save("particle_multiplicity.png");
		
	}
	/*
	public void saveROOT(){
		def ff = new ROOTFile('test0.root');
		ff.addDataSet(hi_efkpfkpfkm_MM_efkpfkpfkm);
		ff.addDataSet(H_FT_e_p_the);
		ff.close();
	}*/

	public void showplots() {

		JFrame frame = new JFrame("SIMULATION");
		frame.setSize(1600, 1000);
		frame.add(myCanvas);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

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
		//System.setProperty("java.awt.headless", "true"); //this needs to uncomment to run script in jlab batch farm
		GStyle.setPalette("kRainBow");
		GStyle.getH1FAttributes().setOptStat("1110");
		//GStyle.getH1FAttributes().setOptStat("000000");
        //GStyle.getFunctionAttributes().setOptStat("1100");

		int count = 0;
		//int maxevents = 1000;
		//int maxevents = 1000000/4;
		long maxevents = 200000000000;
		fdekpkpreqkm_data_v1 ana = new fdekpkpreqkm_data_v1();
		System.out.println(String.format(">>> files from list %s >>>", args[0]));
		String filelist;// = "list_of_files.txt";
		filelist = args[0];
		runType = args[1];
		ArrayList<String> toProcessFileNames = new ArrayList<String>();
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

        //def f1 = new ROOTFile('test.root');
        //def t1 = f1.makeNtuple('h22', 'title', 'x:y:z')

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
        fileWriter1.open("filter_ecpipfpimfp.hipo");
        fileWriter2.open("filter_efpipcpimfp.hipo");

       // */

       	mcoutFile = new File("gen_output_file.txt");
       	mcoutFile.bytes = new byte[0]

       	outFile = new File("output_file.txt");
       	outFile.bytes = new byte[0]

       	// three kaons mixed events for bkg
       	me_outFile = new File("me_output_file.txt");
       	me_outFile.bytes = new byte[0]

       	// two kaons mixed events for bakg
       	me_twoKaon_outFile = new File("me_twoKaon_file.txt");
       	me_twoKaon_outFile.bytes = new byte[0]

       	//public static ArrayList<LorentzVector> lv_my_fdels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
       	lv_my_fdels = new ArrayList<LorentzVector>();
		lv_my_fast_fkps = new ArrayList<LorentzVector>();
		lv_my_slow_fkps = new ArrayList<LorentzVector>();
		lv_my_fkms = new ArrayList<LorentzVector>();

		//public static ArrayList<LorentzVector> lv_my_xi_fdels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms;
		lv_my_xi_fdels = new ArrayList<LorentzVector>();
		lv_my_xi_fast_fkps = new ArrayList<LorentzVector>();
		lv_my_xi_slow_fkps = new ArrayList<LorentzVector>();
		lv_my_xi_fkms = new ArrayList<LorentzVector>();

        for (String runstrg : toProcessFileNames) if(count < maxevents){
        	progresscount++;
            System.out.println(String.format(">>>> file %s >>>>", runstrg));
            File varTmpDir = new File(runstrg);
            if(!varTmpDir.exists()){System.out.println("FILE DOES NOT EXIST");continue;}
            System.out.println("READING NOW " + runstrg);

            if(runstrg == "/Users/akhanal/work/Expt_phys/simulation/xi_1820/rga_spring2019_in_skim14/plot/filter_efkpfkpfkm.hipo"){
            	System.out.println("SPRING 2019 FILE DOES EXIST SETTING Eb = 10.1998 GeV");
            	Eb = 10.1998f;
            	VB = new LorentzVector(0, 0, Eb, Eb);
            } else {
            	Eb = 10.1998f;
            	//Eb = 10.604f;
            	//Eb = 10.604f;
            	VB = new LorentzVector(0, 0, Eb, Eb);
            }

            //Timothy's analysis fitter
            //GenericKinematicFitter research_fitter = new analysis_fitter(Eb);
            //research_fitter = new analysis_fitter(Eb); //research fitter defined as public static GenericKinematicFitter to use Timothys functions
            //PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);

            HipoDataSource reader = new HipoDataSource();
            reader.open(runstrg);
            int filecount = 0;
            while(reader.hasEvent() && count < maxevents) {
            	//DataEvent event = reader.getNextEvent();
            	HipoDataEvent event = reader.getNextEvent();

            	runnum = event.getBank("RUN::config").getInt('run',0);
            	evnum = event.getBank("RUN::config").getInt('event',0);

            	if (runType == "mc"){
					//Perform momentum smearing for MC reconstruction. You don't want to do momentum smearing for real data.
					Random rand = new Random();
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

        //three kaon event mixing
        //public static ArrayList<LorentzVector> lv_my_fdels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms;
        //eventMixing(lv_my_fdels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms);

        //public static ArrayList<LorentzVector> lv_my_xi_fdels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms;
        //twoKaon_eventMixing(lv_my_xi_fdels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps);

        //threeKaon_eventMixing(lv_my_fdels, lv_my_fast_fkps, lv_my_slow_fkps, lv_my_fkms);
        //threeKaon_eventMixing(lv_my_xi_fdels, lv_my_xi_fast_fkps, lv_my_xi_slow_fkps, lv_my_xi_fkms);

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
		//*/		
		System.out.println("Total events : " + count);
		//ana.analyze();
		ana.plot();
		ana.showplots();
		//ana.save();
		//ana.saveROOT();
		
		System.out.println("Good Bye !!!!!");
		println("Ha Ha Ha Good Bye !!!!!");
	}
		
	
	
}
