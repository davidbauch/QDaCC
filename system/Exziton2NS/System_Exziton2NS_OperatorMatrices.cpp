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

        MatrixXcd test1,test2;

        OperatorMatrices() {};
        OperatorMatrices(const Parameters &p) {
            Timer &timer = createTimer("Operator Matrices");
            timer.start();
            generateOperators(p);
            timer.end();
            logs.level2("Done creating operator matrices! Elapsed time is {}ms\n",timer.getWallTime()*1E3);
            outputOperators(p);
        }

        static MatrixXcd tensor(const MatrixXcd &a, const MatrixXcd &b) { //TODO: move to header file
            assert(a.rows() == a.cols() && b.rows() == b.cols());
            MatrixXcd ret = MatrixXcd::Zero(a.cols()*b.cols(), a.rows()*b.rows());
            for (int i = 0; i < a.cols(); i++)
                for (int j = 0; j < a.rows(); j++)
                    for (int k = 0; k < b.cols(); k++) 
                        for (int l = 0; l < b.rows(); l++) {
                            ret(i*a.cols()+k, j*a.rows()+l) = a(i,j)*b(k,l);
                    }
            return ret;
        }

        // Maps the operator op onto N states (e.g. atomic operator onto N photons)
        static MatrixXcd expand_atomic_operator(const MatrixXcd &op, int N) {
            MatrixXcd ret = MatrixXcd::Zero(op.rows()*(N+1),op.rows()*(N+1));
            for (int n = 0; n <= N; n++)
                for (int i = 0; i < op.cols(); i++)
                    for (int j = 0; j < op.rows(); j++)
                        ret(n*op.cols()+i,n*op.rows()+j) = op(i,j);
            return ret;
        }

        // Stretches the operator op onto N states (e.g. photonic operator onto 2 atomic levels)
        static MatrixXcd expand_photonic_operator(const MatrixXcd &op, int N) {
            MatrixXcd ret = MatrixXcd::Zero(op.rows()*N,op.rows()*N);
            for (int i = 0; i < op.cols(); i++)
                for (int j = 0; j < op.rows(); j++)
                    for (int n = 0; n < N; n++)
                        for (int m = 0; m < N; m++)
                            if (n==m)
                                ret(i*N+n,j*N+m) = op(i,j);
            return ret;
        }

        void generateOperators(const Parameters &p) {
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
            //test
            MatrixXcd N1,N2;
            N1 = MatrixXcd::Zero(2,2);
            N2 = MatrixXcd::Zero(p.maxPhotonNumber+1,p.maxPhotonNumber+1);
            N1 <<   0,0,
                    1,0;
            for (int i = 0; i < p.maxPhotonNumber; i++)
                    N2(i+1,i) = sqrt(i+1);
            std::cout << "N1:\n" << N1 << "\nN\n" << N2 << std::endl;
            test1 = expand_atomic_operator(N1,p.maxPhotonNumber);
            test2 = expand_photonic_operator(N2,2);
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
            logs.level2("Hamiltonoperator done! Used:\n{}\nSetting initial rho as pure state with rho_0 = {}... ",H_used,p.rho_0);
            // rho
            rho(p.rho_0,p.rho_0) = 1;
        }

        void outputOperators(const Parameters &p) {
            if (p.outputOperators > 0) {
                std::ostringstream out;
                Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
                if (p.outputOperators > 2) {
                    out << "General Operators:\natom_exited\n" << atom_exited.format(CleanFmt) << std::endl;
                    out << "atom_ground\n" << atom_ground.format(CleanFmt) << std::endl;
                    out << "atom_sigmaplus\n" << atom_sigmaplus.format(CleanFmt) << std::endl;
                    out << "atom_sigmaminus\n" << atom_sigmaminus.format(CleanFmt) << std::endl;
                    out << "atom_inversion\n" << atom_inversion.format(CleanFmt) << std::endl;
                    out << "photon_create\n" << photon_create.format(CleanFmt) << std::endl;
                    out << "photon_annihilate\n" << photon_annihilate.format(CleanFmt) << std::endl;
                    out << "photon_n\n" << photon_n.format(CleanFmt) << std::endl;
                }
                out << "Hamilton and Rho:\nH=H_0+H_I (no RWA)\n" << H.format(CleanFmt) << std::endl;
                out << "H_0\n" << H_0.format(CleanFmt) << std::endl;
                out << "H_I\n" << H_I.format(CleanFmt) << std::endl;
                out << "H_used\n" << H_used.format(CleanFmt) << std::endl;
                out << "rho\n" << rho.format(CleanFmt) << std::endl;
                out << "test1\n" << test1.format(CleanFmt)<< "\ntest2\n" << test2.format(CleanFmt) << std::endl;
                logs.level2(out.str());
                if (p.outputOperators == 3)
                    exit(0);
            }
        }
};