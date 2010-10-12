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



