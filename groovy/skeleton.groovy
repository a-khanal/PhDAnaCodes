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
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import javax.swing.JFrame;

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

    public EmbeddedCanvasTabbed myCanvas;

    public DataGroup dg_kinematics;

    final int BLUE = 9;
    final int LIGHTGREEN = 3;
    final int LIGHTBROWN = 45;
    final int PINK = 46;

    public eppippim(){
    //  NFTElec = 0;
    //  Eb = 10.575f;
        Eb = 10.604f;
    //  Eb = 7.54626f;
    //  Eb = 6.535f
        Mp = (float) PDGDatabase.getParticleById(2212).mass();
    //  Mp = 0.93827f;
        
        VB = new LorentzVector(0, 0, Eb, Eb);
        VT = new LorentzVector(0, 0, 0, Mp);

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
                // Mike's and Alan's momentum correction 
                //float pipCor = FTBmass/(FTBmass + Math.abs(chi2pid) * 0.0271f); //from TOF mass
                //float pipCor = (float) PDGDatabase.getParticleById(211).mass()/(float)(PDGDatabase.getParticleById(211).mass() + chi2pid * 0.0271f);   //from nominal mass
                pip_part_ind = k;
                //pip_mom = mom*pipCor; //with Alan's correction
                //pip_mom = mom;
                //pip_the = the;
                //pip_phi = (float) Math.toDegrees(Math.atan2(py, px));
                //pip_vx = vx;
                //pip_vy = vy;
                //pip_vz = vz;
                //pip_ftb_beta = ftbbe;
                //pip_status = status;
                
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
                                //Vpip = new LorentzVector(pip_px*pipCor, pip_py*pipCor, pip_pz*pipCor, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
                                //Vpip = new LorentzVector(pip_px, pip_py, pip_pz, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
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
                                //Vpipc = new LorentzVector(pip_px*pipCor, pip_py*pipCor, pip_pz*pipCor, Math.sqrt(pip_mom * pip_mom + PDGDatabase.getParticleById(211).mass()*PDGDatabase.getParticleById(211).mass()));
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
                // Mike's and Alan's momentum correction 
                //float pimCor = FTBmass/(FTBmass + Math.abs(chi2pid) * 0.0125f);// from TOF mass
                //float pimCor = (float)PDGDatabase.getParticleById(-211).mass()/(float)( PDGDatabase.getParticleById(-211).mass()+ chi2pid * 0.0125f); // from nominal mass
                pim_part_ind = k;
                //pim_mom = mom;
                //pim_mom = mom*pimCor; //with Alan's correction
                //pim_the = the;
                //pim_phi = (float) Math.toDegrees(Math.atan2(py, px));
                //pim_vx = vx;
                //pim_vy = vy;
                //pim_vz = vz;
                //pim_ftb_beta = ftbbe;
                
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
                //kp_mom = mom;
                //kp_the = the;
                //kp_phi = (float) Math.toDegrees(Math.atan2(py, px));
                //kp_px = px;
                //kp_py = py;
                //kp_pz = pz;
                //kp_vx = vx;
                //kp_vy = vy;
                //kp_vz = vz;
                //kp_ftb_beta = ftbbe;
                
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
                //km_mom = mom;
                //km_the = the;
                //km_phi = (float) Math.toDegrees(Math.atan2(py, px));
                //km_px  = px;
                //km_py  = py;
                //km_pz  = pz;
                //km_vx = vx;
                //km_vy = vy;
                //km_vz = vz;
                //km_ftb_beta = ftbbe;
                
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
                // Mike's and Alan's momentum correction 
                //float proCor1 = FTBmass/(FTBmass + Math.abs(chi2pid) * 0.00893f); //from TOF mass
                //float proCor2 = FTBmass/(FTBmass + Math.abs(chi2pid) * 0.004f);
                //float proCor = proCor1 * proCor2;
                
                //float proCor1 = (float)PDGDatabase.getParticleById(2212).mass()/(float)(PDGDatabase.getParticleById(2212).mass() + chi2pid * 0.00893f); // from nominal mass
                //float proCor2 = (float)PDGDatabase.getParticleById(2212).mass()/(float)(PDGDatabase.getParticleById(2212).mass() + chi2pid * 0.004f);
                //float proCor = proCor1 * proCor2;
                prot_part_ind = k;
                //prot_mom = mom;
                //prot_mom = mom * proCor; //with Alan's correction
                //prot_the = the;
                //prot_phi = (float) Math.toDegrees(Math.atan2(py, px));
                //prot_px = px;
                //prot_py = py;
                //prot_pz = pz;
                //prot_vx = vx;
                //prot_vy = vy;
                //prot_vz = vz;
                //prot_ftb_beta = ftbbe;
                
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

    }
    public void analyze() {

    }

    public void FillHists() {

        if(found_eFD){      
//          H_FT_e_t_f.fill(e_phi, e_the);
//          H_FT_e_p_f.fill(e_phi, e_mom);
//          H_FT_e_p_the.fill(e_the, e_mom);
            H_FT_W_Q2.fill(e_W, e_Q2);
            H_FT_W.fill(e_W);
            H_FT_Q2.fill(e_Q2); 
            H_virphoton.fill(e_virphoton);  
        }

    }

    public void plot(){
        myCanvas = new EmbeddedCanvasTabbed("Kinematics");

        myCanvas.getCanvas("Kinematics").divide(2, 2);
        myCanvas.getCanvas("Kinematics").setSize(1600, 1000);
        myCanvas.getCanvas("Kinematics").setGridX(false);
        myCanvas.getCanvas("Kinematics").setGridY(false);
        myCanvas.getCanvas("Kinematics").setAxisFontSize(18);
        myCanvas.getCanvas("Kinematics").setAxisTitleSize(24);
        myCanvas.getCanvas("Kinematics").draw(dg_kinematics);
        myCanvas.getCanvas("Kinematics").getPad(2).getAxisZ().setLog(true);

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
        int maxevents = 100000;
        //long maxevents = 200000000000;
        eppippim ana = new eppippim();
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
                DataEvent event = reader.getNextEvent();
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
        //ana.analyze();
        ana.plot();
        ana.showplots();
        //ana.plot();
        //ana.save();
        
        System.out.println("Good Bye !!!!!");

    }

 }

