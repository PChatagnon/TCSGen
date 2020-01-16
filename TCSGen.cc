/* 
 * File:   TCSGen.cc
 * Author: rafopar
 *
 * Created on January 11, 2020, 2:53 PM
 */

#include <TF1.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iomanip>
#include <fstream>
#include <TRandom2.h>
#include <TTCSCrs.h>
#include <TTCSKine.h>
#include <KinFunctions.h>
#include <TLorentzVector.h>

#include <cstdlib>
#include <iostream>

using namespace std;
using namespace KinFuncs;

/*
 * 
 */
int main(int argc, char** argv) {

    // ==================================
    // ==== Reading the input config file
    // ==================================

    ifstream inpconfig("GenOptions.dat");

    map<std::string, std::string> m_Settings;
    if (inpconfig.is_open()) {
        while (!inpconfig.eof()) {
            std::string Key;
            std::string Val;
            inpconfig>>Key;
            inpconfig>>Val;
            m_Settings[Key] = Val;
            //cout<<setw(10)<<Key<<setw(20)<<m_Settings[Key]<<endl;
        }
    } else {
        cout << "Can not open the file GenOptions.dat" << endl;
        cout << "So can not initialize settings " << endl;
        cout << "Exiting" << endl;
        exit(1);
    }


    const double PI = 3.14159265358979312;
    const double radian = 57.2957795130823229;
    const double Mp = 0.9383;
    const double Me = 0.00051;


    int Nsim;
    double Eb;
    double t_lim;
    double Eg_min;
    double Eg_max;
    bool isLund;
    double q2_cut;

    // ================== Limits defined by User ==========
    double t_limUser;
    double Eg_minUser;
    double Eg_maxUser;
    double q2_cutUser;
    double Minv_MinUser;
    double Minv_MaxUser;

    // ================== Limits constrained by kimenatics ==========
    // ================== If User limits are out of the kinematic range, then kinematic limits will be used instead
    double t_limKine;
    double Eg_minKine;
    double Eg_maxKine;
    double q2_cutKine;
    double Minv_MinKine;
    double Minv_MaxKine;


    for (map<std::string, std::string>::iterator it = m_Settings.begin(); it != m_Settings.end(); it++) {

        std::string key = (*it).first;
        std::string val = (*it).second;

        if (key.compare("Nsim") == 0) {
            Nsim = atoi(val.c_str());
        } else if (key.compare("Eb") == 0) {
            Eb = atof(val.c_str());
        } else if (key.compare("tLim") == 0) {
            t_lim = atof(val.c_str());
        } else if (key.compare("EgMin") == 0) {
            Eg_minUser = atof(val.c_str());
        } else if (key.compare("EgMax") == 0) {
            Eg_maxUser = atof(val.c_str());
        } else if (key.compare("MinvMax") == 0) {
            Minv_MaxUser = atof(val.c_str());
        } else if (key.compare("MinvMin") == 0) {
            Minv_MinUser = atof(val.c_str());
        } else if (key.compare("q2Cut") == 0) {
            q2_cut = atof(val.c_str());
        } else if (key.compare("LUND") == 0) {
            isLund = atof(val.c_str());
        }

    }

    cout << "Nsim = " << Nsim << endl;
    cout << "Eb = " << Eb << endl;
    cout << "t_lim = " << t_lim << endl;
    cout << "Eg_min = " << Eg_min << endl;
    cout << "Eg_max = " << Eg_max << endl;
    cout << "q2_cut = " << q2_cut << endl;
    cout << "IsLund = " << isLund << endl;


    // =====================================================================
    // ==== We know the beam energy, so Eg_maxKine is Eb
    // =====================================================================

    Eg_maxKine = Eb;
    Eg_max = TMath::Min(Eg_maxUser, Eg_maxKine);

    // ======= Defining Q2, since in most of formulas Q2 will be used instead of Minv
    double Q2MinUser = Minv_MinUser*Minv_MinUser;

    // ======= Now for a givem Q2MinUser, we should check the EgMinUser, if it is below
    // ======= the threshold to produce Min_MinUser, then it should Eg_min should be adjusted
    // ======= to be the threshold for Minv production.
    Eg_minKine = (Q2MinUser + 2*sqrt(Q2MinUser)*Mp)/2*Mp;
    
    Eg_min = TMath::Max(Eg_minKine, Eg_minUser);
    
    cout<<"Eg_min = "<<Eg_min<<endl;
    
    TRandom2 rand;
    rand.SetSeed(0);

    TTCSKine tcs_kin1(Mp, Eb);
    TTCSCrs crs_lmlp;

    TLorentzVector target(0., 0., 0., Mp);
    TLorentzVector Lcm;

    TFile *file_out = new TFile("tcs_gen.root", "Recreate");
    ofstream out_dat("tcs_gen.dat");

    TH2D *h_ph_h_ph_cm1 = new TH2D("h_ph_h_ph_cm1", "", 200, 0., 360., 200, 0., 360.);
    TH2D *h_th_g_th_cm1 = new TH2D("h_th_g_th_cm1", "", 200, 0., 180., 200, 0., 180.);

    //================= Definition of Tree Variables =================
    double Eg, Minv, t, Q2;
    double psf, crs_BH, crs_INT, crs_int;
    double psf_flux, flux_factor;
    TLorentzVector L_em, L_ep, L_prot;
    TLorentzVector L_gprime;

    TTree *tr1 = new TTree("tr1", "TCS MC events");
    tr1->Branch("L_em", "TLorentzVector", &L_em, 3200, 99);
    tr1->Branch("L_ep", "TLorentzVector", &L_ep, 3200, 99);
    tr1->Branch("L_prot", "TLorentzVector", &L_prot, 3200, 99);
    tr1->Branch("Eg", &Eg, "Eg/D");
    tr1->Branch("Q2", &Q2, "Q2/D");
    tr1->Branch("t", &t, "t/D");
    tr1->Branch("psf", &psf, "psf/D");
    tr1->Branch("flux_factor", &flux_factor, "flux_factor/D");
    tr1->Branch("crs_BH", &crs_BH, "crs_BH/D");
    tr1->Branch("crs_INT", &crs_INT, "crs_INT/D");

    for (int i = 0; i < Nsim; i++) {
        if (i % 50000 == 0) {
            cout.flush() << "Processed " << i << " events, approximetely " << double(100. * i / double(Nsim)) << "%\r";
        }

        double psf_Eg = Eg_max - Eg_min;
        Eg = rand.Uniform(Eg_min, Eg_min + psf_Eg);
        flux_factor = N_EPA(Eb, Eg, q2_cut);
        double s = Mp * Mp + 2 * Mp*Eg;
        double t_min = T_min(0., Mp*Mp, Q2MinUser, Mp*Mp, s);
        double t_max = T_max(0., Mp*Mp, Q2MinUser, Mp*Mp, s);
        double psf_t = t_min - TMath::Max(t_max, t_lim);

        //cout<<"t_min = "<<t_min<<"      t_max = "<<t_max<<"    Eg = "<<Eg<<endl;

        if (t_min > t_lim) {
            t = rand.Uniform(t_min - psf_t, t_min);
            double Q2max = 2 * Mp * Eg + t - (Eg / Mp)*(2 * Mp * Mp - t - sqrt(t * t - 4 * Mp * Mp * t)); // Page 182 of my notebook. Derived using "Q2max = s + t - 2Mp**2 + u_max" relation

            double psf_Q2 = Q2max - Q2MinUser;

            Q2 = rand.Uniform(Q2MinUser, Q2MinUser + psf_Q2);

            double u = 2 * Mp * Mp + Q2 - s - t;
            double th_qprime = acos((s * (t - u) - Mp * Mp * (Q2 - Mp * Mp)) / sqrt(Lambda(s, 0, Mp * Mp) * Lambda(s, Q2, Mp * Mp))); //Byukling Kayanti (4.9)
            double th_pprime = PI + th_qprime;

            double Pprime = 0.5 * sqrt(Lambda(s, Q2, Mp * Mp) / s); // Momentum in c.m. it is the same for q_pr and p_pr

            Lcm.SetPxPyPzE(0., 0., Eg, Mp + Eg);
            L_prot.SetPxPyPzE(Pprime * sin(th_pprime), 0., Pprime * cos(th_pprime), sqrt(Pprime * Pprime + Mp * Mp));
            L_gprime.SetPxPyPzE(Pprime * sin(th_qprime), 0., Pprime * cos(th_qprime), sqrt(Pprime * Pprime + Q2));

            double psf_cos_th = 2.; // cos(th):(-1 : 1)
            double psf_phi_cm = 2 * PI;

            double cos_th = rand.Uniform(-1., -1 + psf_cos_th);
            double sin_th = sqrt(1 - cos_th * cos_th);
            double phi_cm = rand.Uniform(0., 0. + psf_phi_cm);

            double El = sqrt(Q2) / 2.; // Energy of lepton in the rest frame of qprime
            double Pl = sqrt(El * El - Me * Me);

            L_em.SetPxPyPzE(Pl * sin_th * cos(phi_cm), Pl * sin_th * sin(phi_cm), Pl*cos_th, El);
            L_ep.SetPxPyPzE(-Pl * sin_th * cos(phi_cm), -Pl * sin_th * sin(phi_cm), -Pl*cos_th, El);

            L_em.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame
            L_ep.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame

            L_em.Boost(L_gprime.BoostVector()); // Move to the CM Frame
            L_ep.Boost(L_gprime.BoostVector()); // Move to the CM Frame

            L_em.Boost(Lcm.BoostVector()); // Move to the Lab Frame
            L_ep.Boost(Lcm.BoostVector()); // Move to the Lab Frame


            L_gprime.Boost(Lcm.BoostVector());
            L_prot.Boost(Lcm.BoostVector());

            double psf_phi_lab = 2 * PI;
            double phi_rot = rand.Uniform(0., psf_phi_lab);

            L_prot.RotateZ(phi_rot);
            L_gprime.RotateZ(phi_rot);
            L_em.RotateZ(phi_rot);
            L_ep.RotateZ(phi_rot);
            tcs_kin1.SetLemLepLp(L_em, L_ep, L_prot);

            h_ph_h_ph_cm1->Fill(phi_cm * TMath::RadToDeg(), tcs_kin1.GetPhi_cm());
            h_th_g_th_cm1->Fill(acos(cos_th) * TMath::RadToDeg(), tcs_kin1.GetTheta_cm());

            psf = psf_t * psf_Q2 * psf_phi_lab * psf_cos_th*psf_phi_cm;

            //crs_lmlp.Set_SQ2t(s, Q2, t);
            crs_BH = crs_lmlp.Eval_BH(s, Q2, t, -1, tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm()); // -1: cros section is not weighted by L/L0

            double eta = Q2 / (2 * (s - Mp * Mp) - Q2);

            // =========== We want to make sure the kinematics is inside the grid of CFFs, otherwise the cross section is not defined
            if (Q2 < 9. && -t < 0.8 && eta < 0.8 && eta > 0.06) {
                crs_INT = crs_lmlp.Eval_INT(s, Q2, t, -1., tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm(), 2.); //the last argumen "1" is the sc_D
                //crs_INT = crs_lmlp.Eval_INT( tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm(), 1.); //the last argumen "1" is the sc_D
            } else {
                crs_INT = 0;
            }

            tr1->Fill();

            //======================== Write LUND file ================================
            double px_em = L_em.Px();
            double py_em = L_em.Py();
            double pz_em = L_em.Pz();
            double px_ep = L_ep.Px();
            double py_ep = L_ep.Py();
            double pz_ep = L_ep.Pz();
            double px_prot = L_prot.Px();
            double py_prot = L_prot.Py();
            double pz_prot = L_prot.Pz();
            //double vz = rand.Uniform(-7.5, 7.5);
            double vz = 0.;

            //============= Write Header ===================
            out_dat << 3 << setw(5) << 1 << setw(5) << 1 << setw(5) << 0 << setw(5) << 0 << setw(15) << crs_BH << setw(15)
                    << crs_INT << setw(15) << 0 << setw(5) << 0 << setw(5) << 0 << endl;
            // =============== WWrite Particles ============
            //====== e- ======
            out_dat << 1 << setw(5) << -1 << setw(5) << 1 << setw(7) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_em << setw(15) << py_em << setw(15)
                    << pz_em << setw(15) << L_em.E() << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(15) << vz << endl;
            //====== e+ ======
            out_dat << 2 << setw(5) << +1 << setw(5) << 1 << setw(7) << -11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_ep << setw(15) << py_ep << setw(15)
                    << pz_ep << setw(15) << L_ep.E() << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(15) << vz << endl;

            //====== proton ======
            out_dat << 3 << setw(5) << +1 << setw(5) << 1 << setw(7) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_prot << setw(15) << py_prot << setw(15)
                    << pz_prot << setw(15) << L_prot.E() << setw(5) << 0 << setw(5) << 0 << setw(5) << 0 << setw(15) << vz << endl;

        } else {
            cout << " |t_min| > |t_lim|" << endl;
            cout << " t_min =  " << t_min << "   t_lim = " << t_lim << "  Eg = " << Eg << endl;
        }
    }

    tr1->Write();
    h_ph_h_ph_cm1->Write();
    h_th_g_th_cm1->Write();


    file_out->Close();

    return 0;
}

