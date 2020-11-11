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
#include <random>

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
	int n_perfile;
	double Eb;
	double t_min;
	double t_max;
	double Eg_min;
	double Eg_max;
	bool isLund;
	double q2_cut;

	// ================== Limits defined by User ==========
	double t_minUser;
	double t_maxUser;
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

	//===================  Target
	double offset_targ;
	double l_targ;
	
	//===================  D-term
	double Dterm;
	
	//=================== Kine cuts
	double Theta_elec_min,Theta_elec_max;
	double P_elec_min;
	double Theta_posi_min,Theta_posi_max;
	double P_posi_min;
	
	//=================== File format
	bool write_lund;
	bool write_root;

	for (map<std::string, std::string>::iterator it = m_Settings.begin(); it != m_Settings.end(); it++) {

		std::string key = (*it).first;
		std::string val = (*it).second;

		if (key.compare("Nsim") == 0) {
			Nsim = atoi(val.c_str());
		} else if (key.compare("Eb") == 0) {
			Eb = atof(val.c_str());
		} else if (key.compare("tMin") == 0) {
			t_minUser = atof(val.c_str());
		} else if (key.compare("tMax") == 0) {
			t_maxUser = atof(val.c_str());
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
		} else if (key.compare("target_l") == 0) {
			l_targ = atof(val.c_str());
		} else if (key.compare("target_off") == 0) {
			offset_targ = atof(val.c_str());
		} else if (key.compare("Dterm") == 0) {
			Dterm = atof(val.c_str());
		} else if (key.compare("n_perfile") == 0) {
			n_perfile = atoi(val.c_str());
		}
		//kine cuts
		 else if (key.compare("Theta_elec_min") == 0) {
			Theta_elec_min = atof(val.c_str());
		} else if (key.compare("Theta_elec_max") == 0) {
			Theta_elec_max = atof(val.c_str());
		} else if (key.compare("P_elec_min") == 0) {
			P_elec_min = atof(val.c_str());
		} else if (key.compare("Theta_posi_min") == 0) {
			Theta_posi_min = atof(val.c_str());
		} else if (key.compare("Theta_posi_max") == 0) {
			Theta_posi_max = atof(val.c_str());
		} else if (key.compare("P_posi_min") == 0) {
			P_posi_min = atof(val.c_str());
		}
		
		// file format
		else if (key.compare("write_lund") == 0) {
			write_lund = (atoi(val.c_str())==1);
		} else if (key.compare("write_root") == 0) {
			write_root = (atoi(val.c_str())==1);
		}

	}

	cout << "Nsim = " << Nsim << endl;
	cout << "Eb = " << Eb << endl;
	cout << "t_min = " << t_minUser << endl;
	cout << "t_max = " << t_maxUser << endl;
	cout << "Eg_min = " << Eg_minUser << endl;
	cout << "Eg_max = " << Eg_maxUser << endl;
	cout << "q2_cut = " << q2_cut << endl;
	cout << "IsLund = " << isLund << endl;
	cout << "target_l = " << l_targ << endl;
	cout << "target_off = " << offset_targ << endl;
	cout << "write_lund = " << write_lund << endl;
	cout << "write_root = " << write_root << endl;

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


	Eg_min =  Eg_minUser;


	TRandom2 rand;
	rand.SetSeed(0);

	TTCSKine tcs_kin1(Mp, Eb);
	TTCSCrs crs_lmlp;

	TLorentzVector target(0., 0., 0., Mp);
	TLorentzVector Lcm;

	TFile *file_out = new TFile("tcs_gen.root", "Recreate");
	ofstream out_dat;
	int file_number=1;
	
	if(write_lund)out_dat.open(Form("tcs_gen_%d.txt", file_number), ofstream::out);("tcs_gen_1.dat");

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

	int i=0;
	
	bool previousEvent = false;	
	while( i < Nsim ){

		if (i % 50000 == 0) {
			cout.flush() << "Processed " << i << " events, approximetely " << double(100. * i / double(Nsim)) << "%\r";
		}

		if( (i)%n_perfile == 0 && previousEvent == true && write_lund){
                             
                                        out_dat.close();
                                        file_number++;
                                        out_dat.open(Form("tcs_gen_%d.txt", file_number), ofstream::out);
                    previousEvent=false;
                }

//cout<<Eg_min<<" "<<Eg_max<<endl;
		double psf_Eg = Eg_max - Eg_min;
		Eg = rand.Uniform(Eg_min, Eg_min + psf_Eg);
		flux_factor = N_EPA(Eb, Eg, q2_cut);
		double s = Mp * Mp + 2 * Mp*Eg;

		bool goodEvent=true;
		if(s<(sqrt(Q2MinUser)+Mp)*(sqrt(Q2MinUser)+Mp)){goodEvent=false;previousEvent=false;continue;}
		
		// cout<<s<<" "<<((sqrt(Q2MinUser)+Mp)*(sqrt(Q2MinUser)+Mp))<<endl;
		double t_min = TMath::Max(T_max(0., Mp*Mp, Q2MinUser, Mp*Mp, s), t_minUser);
		double t_max = TMath::Min(T_min(0., Mp*Mp, Q2MinUser, Mp*Mp, s), t_maxUser);
		double psf_t = t_max - t_min;

		//cout<<"t_min = "<<T_max(0., Mp*Mp, Q2MinUser, Mp*Mp, s)<<"      t_minUser = "<<t_minUser<<endl;
		//cout<<"t_max = "<<T_min(0., Mp*Mp, Q2MinUser, Mp*Mp, s)<<"      t_maxUser = "<<t_maxUser<<endl;
		// cout<<"t_min = "<<t_min<<"      t_max = "<<t_max<<"    Eg = "<<Eg<<" psft "<<psf_t<<endl;

		if (t_min <= t_max && Eg > (Q2MinUser + 2*sqrt(Q2MinUser)*Mp)/2*Mp) {
			t = rand.Uniform(t_min, t_min+psf_t);

			
			//double Q2max = 2 * Mp * Eg + t - (Eg / Mp)*(2 * Mp * Mp - t - sqrt(t * t - 4 * Mp * Mp * t)); // Page 182 of my notebook. Derived using "Q2max = s + t - 2Mp**2 + u_max" relation
			double Q2maxKine = s+Mp*Mp-(1/(2*Mp*Mp))*((s+Mp*Mp)*(2*Mp*Mp-t)-sqrt(Lambda(s,Mp*Mp,0.0)*Lambda(t,Mp*Mp,Mp*Mp)));//equation for s2 p121(5/5.11) in Byckling

			double Q2max = TMath::Min(Q2maxKine, Minv_MaxUser*Minv_MaxUser);
			double psf_Q2 = Q2max - Q2MinUser;

			Q2 = rand.Uniform(Q2MinUser, Q2MinUser + psf_Q2);

			//cout<<Q2<<" "<<psf_Q2<<endl;

			double u = 2 * Mp * Mp + Q2 - s - t;
			double th_qprime = acos((s * (t - u) - Mp * Mp * (Q2 - Mp * Mp)) / sqrt(Lambda(s, 0, Mp * Mp) * Lambda(s, Q2, Mp * Mp))); //Byukling Kayanti (4.9)
			double th_pprime = PI + th_qprime;

			double Pprime = 0.5 * sqrt(Lambda(s, Q2, Mp * Mp) / s); // Momentum in c.m. it is the same for q_pr and p_pr

			Lcm.SetPxPyPzE(0., 0., Eg, Mp + Eg);
			L_prot.SetPxPyPzE(Pprime * sin(th_pprime), 0., Pprime * cos(th_pprime), sqrt(Pprime * Pprime + Mp * Mp));
			L_gprime.SetPxPyPzE(Pprime * sin(th_qprime), 0., Pprime * cos(th_qprime), sqrt(Pprime * Pprime + Q2));

			double psf_cos_th = 2.; // cos(th):(-1 : 1)
			double psf_phi_cm = 2 * PI;

			std::random_device rd;     // only used once to initialise (seed) engine
			std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
			std::uniform_int_distribution<int> uni(0,1); // guaranteed unbiased
			int th= uni(rng); //get 0 or 1
			
			
			
			/*double cos_th;
			if(th)cos_th=rand.Uniform(cos((PI/180.)*50.), cos((PI/180.)*70.));//cos((PI/180.)*60.);
			else cos_th=rand.Uniform(cos((PI/180.)*110.), cos((PI/180.)*130.));*/
			double cos_th = rand.Uniform(-1., -1 + psf_cos_th);
			double sin_th = sqrt(1 - cos_th * cos_th);
			
			/*double phi_cm;
			if(th)phi_cm=rand.Uniform(-50.*(PI/180.), 50.*(PI/180.));
			else {
			phi_cm=rand.Uniform(-50.*(PI/180.), 50.*(PI/180.));
			if(phi_cm<0.)phi_cm=phi_cm+PI;
			else phi_cm=phi_cm-PI;
			}*/
			
			double phi_cm = rand.Uniform(0., 0. + psf_phi_cm);
			//double phi_cm = rand.Uniform(-(1./2.) *PI, PI*(-1./4.));

			//if(!th)cout<<"Phi cm   "<<phi_cm<<endl;

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
			//cout<<crs_BH<<endl;
			double eta = Q2 / (2 * (s - Mp * Mp) - Q2);

			// =========== We want to make sure the kinematics is inside the grid of CFFs, otherwise the cross section is not defined
			if (Q2 < 9. && -t < 0.8 && eta < 0.8 && eta > 0.06) {
				crs_INT = crs_lmlp.Eval_INT(s, Q2, t, -1., tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm(), Dterm); //the last argumen "1" is the sc_D
				//cout<<Dterm<<endl;
				//crs_INT = crs_lmlp.Eval_INT( tcs_kin1.GetPhi_cm(), tcs_kin1.GetTheta_cm(), 1.); //the last argumen "1" is the sc_D
			} else {
				crs_INT = 0;
			}

			

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
			double vz = rand.Uniform(offset_targ-l_targ/2., offset_targ+l_targ/2.);



			if(/*psf<0. ||*/ L_em.Theta()*TMath::RadToDeg()<Theta_elec_min || L_ep.Theta()*TMath::RadToDeg()<Theta_posi_min || L_em.Theta()*TMath::RadToDeg()>Theta_elec_max || L_ep.Theta()*TMath::RadToDeg()>Theta_posi_max || L_em.P()<P_elec_min || L_ep.P()<P_posi_min ){goodEvent=false;previousEvent=false;}

			if(psf<0.0)cout<<"bug in the psf"<<endl;

			double tot_weight = (crs_BH+crs_INT)*flux_factor*psf;
			if(goodEvent){i++;previousEvent=true;};
			if(goodEvent && write_root)tr1->Fill();
			if(goodEvent && write_lund){
				//============= Write Header ===================
				out_dat<<3<<setw(5)<<1<<setw(5)<<1<<setw(15)<<psf<<setw(15)<<crs_BH<<setw(15)<<crs_INT<<setw(15)<<flux_factor<<setw(15)<<crs_INT<<setw(15)<<psf<<setw(15)<<tot_weight<<endl;
				// =============== WWrite Particles ============
				 // Writing Proton
                                        out_dat<<1<<setw(5)<<1<<setw(5)<<1<<setw(7)<<2212<<setw(5)<<0<<setw(5)<<0<<setw(15)<<L_prot.Px()<<setw(15)<<L_prot.Py()<<setw(15)<<L_prot.Pz();
                                        out_dat<<setw(15)<<L_prot.E()<<setw(15)<<Mp<<setw(15)<<0.<<setw(15)<<0.<<setw(15)<<vz<<endl;
                                        // Writing Electron
                                        out_dat<<2<<setw(5)<<-1<<setw(5)<<1<<setw(7)<<11<<setw(5)<<0<<setw(5)<<0<<setw(15)<<L_em.Px()<<setw(15)<<L_em.Py()<<setw(15)<<L_em.Pz();
                                        out_dat<<setw(15)<<L_em.E()<<setw(15)<<Me<<setw(15)<<0.<<setw(15)<<0.<<setw(15)<<vz<<endl;
                                        // Writing Positron
                                        out_dat<<3<<setw(5)<<1<<setw(5)<<1<<setw(7)<<-11<<setw(5)<<0<<setw(5)<<0<<setw(15)<<L_ep.Px()<<setw(15)<<L_ep.Py()<<setw(15)<<L_ep.Pz();
                                        out_dat<<setw(15)<<L_ep.E()<<setw(15)<<Me<<setw(15)<<0.<<setw(15)<<0.<<setw(15)<<vz<<endl;

			}
		} else {
			//cout << " |t_min| > |t_lim|" << endl;
			//cout << " t_min =  " << t_min << "   t_max = " << t_max << "  Eg = " << Eg << endl;
		}
	}

	tr1->Write();
	h_ph_h_ph_cm1->Write();
	h_th_g_th_cm1->Write();


	file_out->Close();

	return 0;
}
