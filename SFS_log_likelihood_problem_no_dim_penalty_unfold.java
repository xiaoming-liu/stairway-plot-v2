/**
 *
 * Copyright (c) @author Xiaoming Liu, Ph.D.
 * Associate Professor,
 * USF Genomics,
 * College of Public Health,
 * University of South Florida at Tampa
 * 
 * This source code is distributed under the Artistic License 2.0
 * 
 * The license can be found at 
 * https://opensource.org/licenses/Artistic-2.0
 */

import java.util.*;
import swarmops.Problem;
import swarmops.Tools;

public class SFS_log_likelihood_problem_no_dim_penalty_unfold extends Problem{
    int n;
    double L;
    double[][] logpara;
    double[] logfac;
    double[] logexi;
    double[] exi;
    double[] xi;
    boolean precomputesReady=false;
    boolean expectedXiReady=false;
    double maxfit;
    double[] lowerBound;
    double[] upperBound;
    String[] parameterName;
    boolean[] obs;
    int nobs;
    boolean setdimpenalty=false;
    double dimpenalty=0;
    boolean setautocorr=false;
    double pautocorr=1;
    int[] groupSplitPoints;
    double[][] groupThetaXiMatrix;
    int[] groupThetaXiMatrixMax;
    double[] xiLogFactorial;

    public SFS_log_likelihood_problem_no_dim_penalty_unfold(int N){
        super();
        n=N;
        logfac= new double[n+1];
        logfac[0]=0;
	for(int i=1;i<=n;i++){
            logfac[i]=logfac[i-1]+Math.log(1.0*i);
        }
        logpara=new double[n][];
        for(int i=1;i<=n-1;i++){
            logpara[i]=new double[n-i+1+1];
        }
        for(int i=1;i<=n-1;i++){
            for(int k=2;k<=n-i+1;k++) logpara[i][k]=logfac[n-i-1]+logfac[n-k]-logfac[n-i-k+1]-logfac[n-1];
        }
        logexi=new double[n];
        exi=new double[n];
        xi=new double[n];
        groupSplitPoints=new int[n];
        for(int i=0;i<n;i++)groupSplitPoints[i]=i+2;//default group: each theta(i) as a group
        obs=new boolean[n];
        for(int i=0;i<n;i++) obs[i]=true;//default observed xi: all xi
        nobs=n;
    }
    public double[] getLogExi(double[] groupThetas){
        if(!expectedXiReady) calculateExpectedXi(groupThetas);
        return logexi;
    }
    public boolean getInitialize2(double[] Theta){
        double totale=0;
        for(int i=1;i<=n-1;i++){
            double exii=0;
            for(int k=2;k<=n-i+1;k++) exii+=Math.exp(logpara[i][k])*Theta[k];
            totale+=exii;
            if (exii==0) exii=Double.MIN_VALUE;
            exi[i]=exii;
            logexi[i]=Math.log(exii);
        }
        if(1-totale>0){
            exi[0]=(1-totale);
            logexi[0]=Math.log(exi[0]);
            expectedXiReady=true;
            return true;
        }
        else {
            
            exi[0]=Double.MIN_VALUE;
            logexi[0]=Math.log(exi[0]);
            return false;
        }
        

    }
    public double preFunction(int n, int i, int j) {
        if (j==2) return 1.0;
        return Math.exp(logfac[n-j+1]+logfac[n-i-1]-logfac[n-i-j+1]-logfac[n-1]);
    }
    public void doPrecomputes() {
        groupThetaXiMatrixMax = new int[n];
        groupThetaXiMatrix = new double[n][];
        for (int i=1;i<=n-1;i++) {
            int jmax = n-i+1;
            int ngroups = groupSplitPoints.length-1;
            int igmax;
            for(igmax=0; igmax<ngroups; igmax++) {
                int j1=groupSplitPoints[igmax+1];
                if (j1 > jmax) break;
            }
            groupThetaXiMatrixMax[i] = igmax;
            groupThetaXiMatrix[i] = new double[igmax+1];
            for(int ig=0; ig<=igmax; ig++) {
                int j0=groupSplitPoints[ig];
                int j1=groupSplitPoints[ig+1];
                double pre0 = preFunction(n,i,j0);
                
                double pre1 = ig==igmax ? 0.0 : preFunction(n,i,j1);
                groupThetaXiMatrix[i][ig] = (pre0-pre1)/i;
            }
        }

        xiLogFactorial = new double[n];
        for(int i=1;i<n;i++) xiLogFactorial[i] = gammln(xi[i]+1.0);

        precomputesReady = true;
    }
    public boolean calculateExpectedXi(double[] groupThetas){
        if (!precomputesReady) doPrecomputes();
        //calculate
        double totale=0;
        for(int i=1;i<=n-1;i++){
            double exii=0;
            for(int ig=0; ig<=groupThetaXiMatrixMax[i]; ig++) {
                exii += groupThetaXiMatrix[i][ig] * groupThetas[ig];
            }
            totale+=exii;
            if (exii==0) exii=Double.MIN_VALUE;
            exi[i] = exii;
            logexi[i]=Math.log(exii);
        }
        if(1-totale>0){
            exi[0]=(1-totale);
            logexi[0]=Math.log(exi[0]);
            expectedXiReady=true;
            return true;
        }
        else {
            
            exi[0]=Double.MIN_VALUE;
            logexi[0]=Math.log(exi[0]);
            return false;
        }
    }
    public void setData(double[] Xi){
        for(int i=0;i<=n-1;i++){
            xi[i]=Xi[i];
        }
        double totall=0;
        L=0;
        for(int i=0;i<=n-1;i++){
        
            L+=xi[i];
        }
        
        double eother=xi[0]/L;
        double oother=xi[0];
        for(int i=1;i<n;i++){
            if(xi[i]>0&&obs[i])totall+=Math.log(xi[i]/L)*xi[i]-gammln(xi[i]+1.0);
            else{
                eother+=xi[i]/L;
                oother+=xi[i];
            }
        }
        totall+=Math.log(eother)*oother-gammln(oother+1.0);
        totall+=gammln(L+1.0);//add constant
        maxfit=totall;//no panelty
        precomputesReady = false;

    }

    public double getLogLikelihood(double[] groupThetas){
        if (!precomputesReady) doPrecomputes();
        calculateExpectedXi(groupThetas);
        double totall=0;
        double eother=exi[0];
        double oother=xi[0];
        for(int i=1;i<n;i++){
            if(obs[i]) totall+=logexi[i]*xi[i]-xiLogFactorial[i];
            else{
                eother+=exi[i];
                oother+=xi[i];
            }
        }
        totall+=Math.log(eother)*oother-gammln(oother+1.0);
        totall+=gammln(L+1.0);//add constant

        double adjustlogL=totall;//no panelty
        
        if(setdimpenalty)adjustlogL=totall-dimpenalty;
        return adjustlogL;
    }
    /**
     *
     * @param eXi, expected freq of \xi, from \xi_0 to \xi_n-1
     * @param Xi, observed counts of \xi from \xi_0 to \xi_n-1
     * @return logL
     */
    public double getLogLikelihood(double[] eXi,double[] Xi){
        n=eXi.length;
        for(int i=0;i<=n-1;i++)logexi[i]=Math.log(eXi[i]);
        L=0;
        for(int i=0;i<=n-1;i++){
            xi[i]=Xi[i];
            L+=xi[i];
        }
        double totall=0;
        double eother=Math.exp(logexi[0]);
        double oother=xi[0];
        for(int i=1;i<n;i++){
            if(obs[i]) totall+=logexi[i]*xi[i]-gammln(xi[i]+1.0);
            else{
                eother+=Math.exp(logexi[i]);
                oother+=xi[i];
            }
        }
        totall+=Math.log(eother)*oother-gammln(oother+1.0);
        totall+=gammln(L+1.0);//add constant

        double adjustlogL=totall;//no panelty
        
        if(setdimpenalty)adjustlogL=totall-dimpenalty;
        return adjustlogL;
    }
    public void setDimPenalty(double p){
        dimpenalty=p;
        setdimpenalty=true;
    }
    public void setAutoCorr(double p){
        setautocorr=true;
        pautocorr=p;
    }
    public void setThetaGroup(int[][] Group){
        groupSplitPoints = new int[Group.length+1];
        for(int i=0;i<Group.length;i++)groupSplitPoints[i]=Group[i][0];
        int ilast=Group.length-1;
        int jlast=Group[ilast].length-1;
        groupSplitPoints[ilast+1]=Group[ilast][jlast]+1;
        precomputesReady = false;
    }
    /**
     * set whether xi(i) is used for analysis
     * @param Obs
     */
    public void setObsXi(boolean[] Obs){
        obs=(boolean[])Obs.clone();
        nobs=1;
        for(int i=1;i<n;i++) if(obs[i])nobs++;
        precomputesReady = false;
    }
        public String getName() {
		return "SFS_log_likelihood_problem";
	}
        public int getDimensionality() {
		return groupSplitPoints.length-1;
	}

        public double[] getLowerBound() {
            lowerBound = new double[getDimensionality()];
            Arrays.fill(lowerBound,0);
            
		return lowerBound;
	}

        public double[] getUpperBound() {
            upperBound = new double[getDimensionality()];
            Arrays.fill(upperBound,0.2);
            return upperBound;
	}
        @Override
        public double[] getLowerInit() {
		return getLowerBound();
	}

	@Override
	public double[] getUpperInit() {
		return getUpperBound();
	}
        public double getMinFitness() {
            return -maxfit;
	}
        @Override
	public double getAcceptableFitness() {

            return -maxfit;
	}


        @Override
	public String[] getParameterName() {
            parameterName =new String[getDimensionality()];
            for(int i=0;i<n-1;i++) parameterName[i]="theta"+(i+2);
		return parameterName;
	}
        @Override
	public double fitness(double[] x) {
		assert x != null && x.length == getDimensionality();
        double[] groupThetas = x;
        return -getLogLikelihood(groupThetas);
	}
        @Override
	public boolean enforceConstraints(double[] x) {
		// Enforce boundaries.
		Tools.bound(x, getLowerBound(), getUpperBound());

		// Return feasibility.
		return isFeasible(x);
	}
        @Override
	public boolean isFeasible(double[] x) {//isFeasible come first before geting fitness
		assert x != null && x.length == getDimensionality();
                boolean feasible=true;

                double[] groupThetas = x;
                for (int ig=0; ig<groupThetas.length; ig++) if (groupThetas[ig]==0) return false;

                if(setautocorr){
                    double p2=(1-pautocorr)/2;
                    //assume theta' drawn from an exponential distribution with mean theta, from Bayesian skyline plot
                    for (int ig=1; ig<groupThetas.length-0; ig++) {
                        double lambda=1/groupThetas[ig-1];
                        double cdf=1-Math.exp(-lambda*groupThetas[ig]);
                        if(cdf<p2||cdf>1-p2) return false;
                    }
                    for (int ig=0; ig<groupThetas.length-1; ig++) {
                        double lambda=1/groupThetas[ig+1];
                        double cdf=1-Math.exp(-lambda*groupThetas[ig]);
                        if(cdf<p2||cdf>1-p2) return false;
                    }
                }

                double k=0;
                for (int ig=0; ig<groupThetas.length; ig++) {
                    int i0 = groupSplitPoints[ig];
                    int i1 = groupSplitPoints[ig+1];
                    for (int i=i0; i<i1; i++) k += groupThetas[ig]/(i-1);
                }
                if(k>1)feasible=false;

		return feasible;
	}
        /**
       * Compute LnGamma(x)
       */
	private static double gammln(double xx)
	{
		int j;
		double temp;
		double cof[] = new double[7];
		double stp, half, one, fpf, x, tmp, ser;
		cof[1] = 76.18009172947146;
		cof[2] = -86.50532032941677;
		cof[3] = 24.01409824083091;
        cof[4] = -1.231739572450155;
        cof[5] = 0.001208650973866179;
        cof[6] = -0.000005395239384953;
        stp = 2.5066282746310005;
        half = 0.5;
        one = 1.0;
        fpf = 5.5;
        x = xx ;
        tmp = x + fpf;
        tmp = (x + half) * Math.log(tmp) - tmp;
        ser = 1.000000000190015;
        for (j = 1; j <= 6; j++)
		{
            x = x + one;
            ser = ser + cof[j] / x;
        }
        temp = tmp + Math.log(stp * ser/xx);
		return temp;
	}
    public static void main(String[] args){
        
    }
}
