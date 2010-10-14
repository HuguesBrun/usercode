#include "functions.h"


int* findLeadAndTrail(TLorentzVector* P, int taille){
	float Et = 0;
	float EtMax = 0;
	float EtSec = 0;
	int iteMax = -1;
	int iteSec = -1;
	int * leResu = new int[2];
	for (int i = 0 ; i < taille ; i++){
		Et = P[i].Et();
		if (Et > EtMax) {
			EtMax = Et;
			iteMax = i;
		}
        }
	for (int i = 0; i < taille ; i++){
		Et = P[i].Et();
		if ((Et > EtSec)&&(i!=iteMax)){
			EtSec = Et;
			iteSec = i;
		}
	}
	leResu[0] = iteMax;
	leResu[1] = iteSec;
	return leResu;
}

double CosThetaStar(TLorentzVector p1, TLorentzVector p2){
	TLorentzVector p = p1 + p2;
	TVector3 theBoost = p.BoostVector();
	TVector3 bostDir;
	if ( theBoost.Mag() != 0 ) bostDir = theBoost.Unit(); // / theBoost.Mag());
	else return -1;
	p1.Boost(-theBoost);
	if (p1.Vect().Mag()!=0) return p1.Vect().Dot(bostDir) / p1.Vect().Mag();
	else return -1;	
}

