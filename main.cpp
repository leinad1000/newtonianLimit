
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

// Åpner objekt som skriver til fil
ofstream ofile;

int main(int argc, char *argv[])
{

    int n, N;
    char *outfilename;

    // Reading in outputfile. Abort if there are too few command line arguments

    if (argc<3){
        cout << "Bad usage: " << argv[0] << "read also outputfile and 'n' on the same line"<<endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
        n = atoi(argv[2]); // Ikke større enn 1000, hvis ikke tar det for lang tid
    }

    // Åpner fil
    ofile.open(outfilename);

    // Konstanter
    double omegaNull, omegaMax, omegaStep, deltaNull, deltaMax, deltaDotMax, deltaDotDotMax, thetaNull, phiNull, thetaDotDot, phiDot, phiDotPlussEn, xiMin, xiMax, step, treHalve, femHalve, tjuefemHalve, thetaFemHalve, stepStep, stepHalve, toStep, stepStepHalve;
    // Velger intervall for egenverdier
    omegaNull = 0;
    omegaMax = 6; // Jo høyere denne er, dess flere egenverdier finner du
    N = 1E3; // Høy N først for å sørge for at jeg finner alle egenverdiene i det valgte intervallet.
             // Denne senkes siden når egenverdiene skal finnes med høyere presisjon.
             // Den senkes for å speede opp koden.
    omegaStep = (omegaMax - omegaNull)/((double) (N-1));
    treHalve = 3.0/2.0;
    femHalve = 5.0/2.0;
    tjuefemHalve = 25.0/2.0;
    thetaNull = 1.0;
    phiNull = 0;
    deltaNull = 0;
    deltaMax = 1.0;
    xiMax = 3.6535981;
    xiMin = 1e-6;
    step = (xiMax - xiMin)/((double) (n-1));
    stepStep = pow(step,2.0);
    stepHalve = step/2.0;
    toStep = 2.0*step;
    stepStepHalve = (pow(step,2.0)/2.0);


    // Vektorer
    double *theta = new double[n];
    double *thetaDot = new double[n];
    double *phi = new double[n];
    double *xi = new double[n];
    double *k = new double[n];
    double *l = new double[n];
    double *m = new double[n];
    double *delta = new double[n];
    double *deltaNullVektor = new double[N];
    double *omega = new double[N];

    // Initialverdier
    xi[0] = xiMin;
    theta[0] = thetaNull;
    phi[0] = phiNull;
    phiDot = pow(theta[0],treHalve)*pow(xi[0],2.0);

    // løser Lane-Emden
    for (int i=0;i<n-1;i++){
        xi[i+1] = step*(i+1);
        thetaDot[i] = -phi[i]/(pow(xi[i],2.0));
        thetaDotDot = -(2.0*thetaDot[i])/xi[i] - pow(theta[i],treHalve);
        theta[i+1] = theta[i] + step*thetaDot[i] + stepStepHalve*thetaDotDot;
        phiDotPlussEn = pow(theta[i+1],treHalve)*pow(xi[i+1],2.0);
        phi[i+1] = phi[i] + stepHalve*phiDot + stepHalve*phiDotPlussEn;
        phiDot = phiDotPlussEn;
    }
    thetaDot[n-1] = -phi[n-1]/(pow(xi[n-1],2.0));

    // Løser den radielle egenverdiligningen

    // Finner k og l
    for (int i=0;i<n;i++){
        thetaFemHalve = pow(theta[i],femHalve);
        k[i] = 5.0*thetaFemHalve;
        l[i] = (10.0*thetaFemHalve)/xi[i] + tjuefemHalve*pow(theta[i],treHalve)*thetaDot[i];
    }

    // Finner omega'ene
    for (int i=0;i<N;i++){
        omega[i] = omegaNull + omegaStep*i;
    } 

    // Finner de forkjellige delta-null'ene
    for (int j=0;j<N;j++){
        for (int i=0;i<n;i++){
            m[i] = pow(omega[j],2.0)*pow(theta[i],treHalve) - (5.0*pow(theta[i],treHalve)*thetaDot[i])/xi[i] - (10.0*pow(theta[i],femHalve))/(pow(xi[i],2.0));
        }

        delta[n-1] = deltaMax; // Første initialverdi
        deltaDotMax = (5*thetaDot[n-1]/xi[n-1] - pow(omega[j],2.0))/(tjuefemHalve*thetaDot[n-1]);
        deltaDotDotMax = -(l[n-1]*deltaDotMax + m[n-1]*deltaMax)/k[n-1];
        delta[n-2] = deltaMax - step*deltaDotMax + stepStepHalve*deltaDotDotMax; // Andre initialverdi

        for (int i=n-1;i>1;i--){
            double A = 4 - 2*stepStep*(m[i]/k[i]);
            double B = 2 - step*(l[i]/k[i]);
            double C = 2 + step*(l[i]/k[i]);
            double D = 2 - step*(l[i]/k[i]);
            delta[i-2] = (A/B)*delta[i-1] - (C/D)*delta[i];
        }
        deltaNullVektor[j] = delta[0];
    }

    // Konstanter
    double fraction, forrigeDelta, Delta;
    int antallEgenverdier, teller;
    fraction = 1.0;
    forrigeDelta = deltaNullVektor[0];
    antallEgenverdier = 0;

    // Teller antall egenverdier
    for (int i=0;i<N;i++) {
        Delta = deltaNullVektor[i];
        fraction = Delta/forrigeDelta;
        if (fraction < 0) {
            antallEgenverdier += 1;
        }
        forrigeDelta = Delta;
    }

    // Lager matrise som skal inneholde de forskjellige OmegaMax og OmegaMin omkring hver egenverdi
    double **omegaMatrise = new double*[antallEgenverdier];
    for(int i = 0; i < antallEgenverdier; i++){
        omegaMatrise[i] = new double[2];
    }

    // Fyller opp omegaMatrise
    forrigeDelta = deltaNullVektor[0];
    teller = 0;

    for (int i=1;i<N;i++) {
        Delta = deltaNullVektor[i];
        fraction = Delta/forrigeDelta;
        if (fraction < 0) {
            omegaMatrise[teller][0] = omega[i-1];
            omegaMatrise[teller][1] = omega[i];
            teller += 1;
        }
        forrigeDelta = Delta;
    }

    // Konstanter
    double epsilon, besteDeltaNull;
    epsilon = 1E-3;
    double *egenverdiVektor = new double[antallEgenverdier];
    double **deltaMatrise = new double*[antallEgenverdier];
    for(int i = 0; i < antallEgenverdier; i++){
        deltaMatrise[i] = new double[n];
    }

    N = 1E1; // Senker N for å speede opp koden. Presisjonen blir ikke dårligere, siden det er epsilon som bestemmer presisjonen.
             // Lavere N gjør bare at while-løkken nedenfor entres flere ganger. Dette gjør ikke noe, siden for ved høy N,
             // entres den allerede få ganger. Lavere N gjør at den bare entres maks et par ganger ekstra.

    // Finner egenverdiene innenfor ønsket presisjon
    for (int i=0;i<antallEgenverdier;i++) {
        omegaNull = omegaMatrise[i][0];
        omegaMax = omegaMatrise[i][1];
        omegaStep = (omegaMax - omegaNull)/((double) (N-1));
        besteDeltaNull = 2*epsilon;
        while (fabs(besteDeltaNull) > epsilon) {
            for (int j=0;j<N;j++){
                omega[j] = omegaNull + omegaStep*j;
            }
            for (int j=0;j<N;j++){
                for (int ii=0;ii<n;ii++){
                    m[ii] = pow(omega[j],2.0)*pow(theta[ii],treHalve) - (5.0*pow(theta[ii],treHalve)*thetaDot[ii])/xi[ii] - (10.0*pow(theta[ii],femHalve))/(pow(xi[ii],2.0));
                }

                delta[n-1] = deltaMax; // Første initialverdi
                deltaDotMax = (5*thetaDot[n-1]/xi[n-1] - pow(omega[j],2.0))/(tjuefemHalve*thetaDot[n-1]);
                deltaDotDotMax = -(l[n-1]*deltaDotMax + m[n-1]*deltaMax)/k[n-1];
                delta[n-2] = deltaMax - step*deltaDotMax + stepStepHalve*deltaDotDotMax; // Andre initialverdi

                for (int q=n-1;q>1;q--){
                    double A = 4 - 2*stepStep*(m[q]/k[q]);
                    double B = 2 - step*(l[q]/k[q]);
                    double C = 2 + step*(l[q]/k[q]);
                    double D = 2 - step*(l[q]/k[q]);
                    delta[q-2] = (A/B)*delta[q-1] - (C/D)*delta[q];
                }
                deltaNullVektor[j] = delta[0];
            }
            forrigeDelta = deltaNullVektor[0];
            for (int j=1;j<N;j++) {
                Delta = deltaNullVektor[j];
                fraction = Delta/forrigeDelta;
                if (fraction < 0) {
                    besteDeltaNull = forrigeDelta;
                    omegaNull = omega[j-1];
                    omegaMax = omega[j];
                    omegaStep = (omegaMax - omegaNull)/((double) (N-1));
                }
                forrigeDelta = Delta;
            }

        }
        egenverdiVektor[i] = omegaNull;
        //cout << pow(egenverdiVektor[i],2.0) << endl; // Skriver ut egenverdiene
    }

    cout << "Antall Egenverdier = " << antallEgenverdier << endl;

    // Skriver xi til fil
    for (int i=n-1;i>0;i--) {
        ofile << setw(15) << setprecision(8) << xi[i]/xi[n-1] << "\t";
    }
    ofile << setw(15) << setprecision(8) << xi[0]/xi[n-1] << endl;

    // Regner ut egenvektorene og skriver de til fil
    double Omega;

    for (int i=0;i<antallEgenverdier;i++) {
        Omega = egenverdiVektor[i];
        for (int j=0;j<n;j++){
            m[j] = pow(Omega,2.0)*pow(theta[j],treHalve) - (5.0*pow(theta[j],treHalve)*thetaDot[j])/xi[j] - (10.0*pow(theta[j],femHalve))/(pow(xi[j],2.0));
        }
        delta[n-1] = deltaMax; // Første initialverdi
        deltaDotMax = (5*thetaDot[n-1]/xi[n-1] - pow(Omega,2.0))/(tjuefemHalve*thetaDot[n-1]);
        deltaDotDotMax = -(l[n-1]*deltaDotMax + m[n-1]*deltaMax)/k[n-1];
        delta[n-2] = deltaMax - step*deltaDotMax + stepStepHalve*deltaDotDotMax; // Andre initialverdi

        // Skriver siste verdi til fil, først
        ofile << setw(15) << setprecision(8) << delta[n-1] << "\t";

        for (int j=n-1;j>1;j--){
            double A = 4 - 2*stepStep*(m[j]/k[j]);
            double B = 2 - step*(l[j]/k[j]);
            double C = 2 + step*(l[j]/k[j]);
            double D = 2 - step*(l[j]/k[j]);
            delta[j-2] = (A/B)*delta[j-1] - (C/D)*delta[j];
            ofile << setw(15) << setprecision(8) << delta[j-1] << "\t";
        }
        // Skriver første verdi til fil, sist
        ofile << setw(15) << setprecision(8) << delta[0] << endl;

    }

    // Stenger fil
    ofile.close();

    // Sletter vektorer
    delete [] theta;
    delete [] thetaDot;
    delete [] phi;
    delete [] xi;
    delete [] k;
    delete [] l;
    delete [] m;
    delete [] delta;
    delete [] deltaNullVektor;
    delete [] omega;
    delete [] egenverdiVektor;

    // Sletter matrise
    for(int i = 0; i < 2; ++i) {
        delete [] omegaMatrise[i];
    }

    delete [] omegaMatrise;

    return 0;
}
