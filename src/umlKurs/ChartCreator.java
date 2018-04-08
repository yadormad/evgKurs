package umlKurs;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import static java.lang.Math.*;

public class ChartCreator {
    private double length, D, time;
    private int tn, xn, nSum;
    private double[] lamdas;
    private static final double H = 0.005;

    public ChartCreator(double length, double d, double time, int N, int xn, int tn) {
        Scanner scanner = null;
        try {
            scanner = new Scanner(new File("lamdas.txt"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        this.length = length;
        D = d;
        this.time = time;
        this.tn = tn;
        this.xn = xn;
        this.nSum = N;
        lamdas = new double[nSum];
        for(int i = 0; i < nSum; i++){
            //String string = scanner.nextLine();
            assert scanner != null;
            lamdas[i] = Double.parseDouble(scanner.nextLine());
        }
    }


    public double getLength() {
        return length;
    }

    public double getTime() {
        return time;
    }

    public int getTn() {
        return tn;
    }

    public int getXn() {
        return xn;
    }

    public int getnSum() {
        return nSum;
    }

    public double ant(int i) {
        return ((1 / lamdas[i]) * Math.sin((length * lamdas[i]) / 2)) / ((length / 2)+((length / (4 * lamdas[i])) * Math.sin(2 * length * lamdas[i])));

        /*double alphCPiNL2;
        double divKDC = gammaK - D/c;
        alphCPiNL2 = (alpha/c)*pow(PI*i/length, 2.0);
        double an = divKDC*exp(-gammaK*t)/(alphCPiNL2 - divKDC);
        an += -divKDC*exp(-(alphCPiNL2 + D/c)*t)/(alphCPiNL2 - divKDC);
        an += -exp(-(alphCPiNL2 + D/c)*t);
        return an;*/
    }



    public double decompositionOfOne(int i, double x, double t) { //базис
        double fx = Math.exp((-D)*lamdas[i]*lamdas[i]*t);
        double cosBlaBla = cos(lamdas[i]*x);
        return fx*cosBlaBla;
    }

    public double getSum(double x, double t) {
        double result = 0;
        double tmp1, tmp2;
        for (int i = 1; i <= nSum; i++) {
            tmp1 = ant(i-1);
            tmp2 = decompositionOfOne(i-1, x, t);
            //result += (ant(i) /*+ exp(-gammaK*t)*/)*decompositionOfOne(i, x, t) ;
            result += tmp1*tmp2;
        }
        return result;
    }

    public double[][] composeSeparChartArray(){
        int xnSepar = 300;
        double[][] result = new double[xnSepar +1][tn +1];
        for(int i = 0; i <= xnSepar; i++) {
            for (int j = 0; j <= tn; j++) {
                result[i][j] = getSum(length*i/ xnSepar, time*j/ tn) + 0.5;
            }
        }
        /*for (int j = 0; j <= tn; j++) {
            result[0][j] = exp(-gammaK*time*j/tn);
            result[xn][j] = result[0][j];
        }*/
        return result;
    }

    /*private int calcN(double eps) {
        int startN = 50;
        int step = 10;
        double ans1, ans2;
        do {
            ans1 = 0;
            for(int i = 1; i <= startN; i++) {
                ans1 += ant(i, 0) * decompositionOfOne(i, length/2*startN);
            }
            ans2 = 0;
            for(int i = 1; i <= startN + step; i++) {
                ans1 += ant(i, 0) * decompositionOfOne(i, length/2*(startN + step));
            }
            startN += step;
        } while (abs(ans1 - ans2) > eps);
        return startN;
    }*/

    /*private int calcNx(double errT){
        int startN = 2;
        ArrayList<Double> ans1;// = new ArrayList<>();
        double vx0;
        double step;
        do {
            step = length/startN;
            ans1 = new ArrayList<>();
            for (int j = 0; j <= startN; j++) {
                vx0 = 0;
                for (int i = 1; i <= nSum; i++) {
                    vx0 += ant(i, 0) * decompositionOfOne(i, startN*step);
                }
                ans1.add(vx0);
            }

            for (int j = 0; j < startN; j++) {
                vx0 = 0;
                for (int i = 1; i <= nSum; i++) {
                    vx0 += ant(i, 0) * decompositionOfOne(i, startN*step);
                }
                ans1.add(vx0);
            }

            ans2 = 0;
            for(int i = 1; i <= startN + step; i++) {
                ans1 += ant(i, 0) * decompositionOfOne(i, length/2*(startN + step));
            }
        } while (abs(ans1 - ans2) > errT);
        return startN;
    }*/

    private double[] getNextKStep(double[] q, int k, double gam, double r, double p, double p0, double r0, double pI) {
        double[] nextCol = new double[q.length];
        double[] alf = new double[xn + 1];
        double[] bet = new double[xn + 1];
        alf[1] = -r0/p0;
        bet[1] = q[0]/p0;
        for (int j = 2; j <= xn; j++) {
            alf[j] = -r/(p + r*alf[j-1]);
            bet[j] = (q[j-1] - r*bet[j-1])/(p + r*alf[j-1]);
        }
        nextCol[xn] = (q[xn] - r*bet[xn])/(pI + r0*alf[xn]);
        for (int m = xn - 1; m >= 0; m--) {
            nextCol[m] = alf[m+1]*nextCol[m+1] + bet[m+1];
        }
        return nextCol;
    }

    public double[][] composeRunChartArray(){
        double[][] result = new double[xn +1][tn +1];
        double[] q = new double [xn + 1];
        double gam = (time/ tn)/(Math.pow(length/ xn, 2.0));
        double r = -D*gam;
        double p = 2*gam*D + 1;
        double p0 = 2*gam*D + 1;
        double r0 = 2*r;
        double  pI = p + 2*H*length/ xn *(-r);
        for(int i = 0; i <= xn /2; i++) {
            result[i][0] = 1;
            q[i] = 1;
        }
        for(int i = xn /2 + 1; i <= xn; i++) {
            result[i][0] = 0;
            q[i] = 0;
        }
        for(int k = 0; k < tn; k++) {
            q = getNextKStep(q, k, gam, r, p, p0, r0, pI);
            for(int i = 0; i <= xn; i++) {
                result[i][k+1] = q[i];
            }
        }
        return result;
    }
}
