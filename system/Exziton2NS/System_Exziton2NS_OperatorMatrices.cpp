class OperatorMatrices {
    public:
        MatrixXcd H;
        MatrixXcd H_0;
        MatrixXcd H_I;
        MatrixXcd rho;
        MatrixXcd H_used;
        MatrixXcd photon_create;
        MatrixXcd photon_annihilate;
        MatrixXcd atom_exited;
        MatrixXcd atom_ground;
        MatrixXcd photon_n;
        MatrixXcd atom_sigmaplus;
        MatrixXcd atom_sigmaminus;
        MatrixXcd atom_inversion;

        OperatorMatrices() {};
        OperatorMatrices(const Parameters &p) {
            Timer &timer = createTimer("Operator Matrices");
            timer.start();

            logs.level2("Creating operator matrices, dimension = {}\nCreating base matrices... ",p.maxStates);
            H = MatrixXcd::Zero(p.maxStates,p.maxStates);
            H_0 = MatrixXcd::Zero(p.maxStates,p.maxStates);
            H_I = MatrixXcd::Zero(p.maxStates,p.maxStates);
            photon_create = MatrixXcd::Zero(p.maxStates,p.maxStates);
            photon_annihilate = MatrixXcd::Zero(p.maxStates,p.maxStates);
            photon_n = MatrixXcd::Zero(p.maxStates,p.maxStates);
            atom_exited = MatrixXcd::Zero(p.maxStates,p.maxStates);
            atom_ground = MatrixXcd::Zero(p.maxStates,p.maxStates);
            atom_sigmaplus = MatrixXcd::Zero(p.maxStates,p.maxStates);
            atom_sigmaminus = MatrixXcd::Zero(p.maxStates,p.maxStates);
            atom_inversion = MatrixXcd::Zero(p.maxStates,p.maxStates);
            H_used = MatrixXcd::Zero(p.maxStates,p.maxStates);
            rho = MatrixXcd::Zero(p.maxStates,p.maxStates);

            logs.level2("Done! Finalizing matrices... ");
            for (int n = 0; n < p.maxStates; n++)
                for (int m = 0; m < p.maxStates; m++) {
                    if (n==m) photon_n(n,m) = (int)((n)/2);                       
                        photon_create(n,m) = std::sqrt((int)(m/2+1))*delta(n%2,m%2)*delta((int)(n/2),(int)(m/2)+1);          
                        photon_annihilate(n,m) = std::sqrt((int)(m/2))*delta(n%2,m%2)*delta((int)(n/2),(int)(m/2)-1);
                        atom_sigmaplus(n,m) = delta((int)(n/2),(int)(m/2))*delta((int)(n%2),1)*delta((int)(m%2),0);
                        atom_sigmaminus(n,m) = delta((int)(n/2),(int)(m/2))*delta((int)(n%2),0)*delta((int)(m%2),1);
                        atom_ground(n,m) = delta((int)(n/2),(int)(m/2))*delta((int)(n%2),0)*delta((int)(m%2),0);
                        atom_exited(n,m) = delta((int)(n/2),(int)(m/2))*delta((int)(n%2),1)*delta((int)(m%2),1);
                    }
                    for (int n = 0; n < p.maxStates; n++) {
                        atom_inversion(n,n) = (n%2 == 0) ? -1 : 1;
                    }
            logs.level2("Done! Creating Hamiltonoperator... ");
            /* All possible Hamiltonions */
            // H_0
            H_0 = p.omegaE/2.0*atom_inversion + p.omegaC*photon_n;
            // H_I
            // RWA
            if (p.doRWA) {
                logs.level2("using RWA... ");
                H_I = p.g*( atom_sigmaplus*photon_annihilate + atom_sigmaminus*photon_create );
            }
            // non RWA
            if (!p.doRWA) {
                logs.level2("NOT using RWA... ");
                //H_I = p.g*( (atom_sigmaplus+atom_sigmaminus)*(photon_annihilate+photon_create) );
                H_I = p.g*( atom_sigmaplus*photon_create + atom_sigmaplus*photon_annihilate + atom_sigmaminus*photon_create + atom_sigmaminus*photon_annihilate );
            }
            // H
            H = H_0 + H_I;
            
            if (p.doInteractionPicture) {
                logs.level2("using interaction picture... ");
                H_used = H_I;
            }
            
            if (!p.doInteractionPicture) {
                logs.level2("NOT using interaction picture... ");
                H_used = H;
            }
            logs.level2("Hamiltonoperator done!\nSetting initial rho as pure state with rho_0 = {}... ",p.rho_0);
            // rho
            rho(p.rho_0,p.rho_0) = 1;
            timer.end();
            logs.level2("Done creating operator matrices! Elapsed time is {}ms\n",timer.getWallTime()*1E3);
        }
};