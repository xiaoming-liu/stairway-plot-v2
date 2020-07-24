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
import java.io.*;

import swarmops.*;
import swarmops.optimizers.*;

public class Stairway_unfold_training_testing7 {
     public static void main(String[] args)throws Exception{
        
        boolean verbose=false;
        boolean fitnesstrail=true;
        int numRuns = 1;
        int totalRuns=0;
        int dimFactor = 2000;
        double pautocorr=0.99;    
        if(args.length!=4&&args.length!=3){
            System.out.println("Usage: java Stairway_unfold_training_testing7 <input file name> <numBreaks> <percentage_training> [random_seed]");
            System.exit(0);
        }
        
        
        String infile =args[0];
        int numBreaks=Integer.parseInt(args[1]);
        
        double ptraining=Double.parseDouble(args[2]);
        long seed=(new Random()).nextLong();
        if(args.length==4)seed=Long.parseLong(args[3]);
        PrintWriter out=new PrintWriter(new FileWriter(infile+"."+numBreaks+"_"+ptraining+".addTheta"),true);
        BufferedReader in=new BufferedReader(new FileReader(infile));
        int nline=0;
        while(in.ready()){
            in.readLine();
            nline++;
        }
        in.close();
        in=new BufferedReader(new FileReader(infile));
        int nsfs=nline/2;
        for(int ii=0;ii<nsfs;ii++){
            String line=in.readLine();
            
            StringTokenizer t=new StringTokenizer(line);
            int nt=t.countTokens();
            String[] temp=new String[nt];
            for(int i=0;i<nt;i++)temp[i]=t.nextToken();
            String popid=temp[0];
            int nseq=Integer.parseInt(temp[1]);
            
            double L=Double.parseDouble(temp[2]);
            if(nt>5)L=Double.parseDouble(temp[3]);
            int xibegin=Integer.parseInt(temp[3]);//begin of observed Xi, for example 2 for xi(2)
            if(nt>5)xibegin=Integer.parseInt(temp[4]);//begin of observed Xi, for example 2 for xi(2)
            int xiend=Integer.parseInt(temp[4]);//end of observed Xi, for example n-2 for xi(n-2)
            if(nt>5)xiend=Integer.parseInt(temp[5]);//end of observed Xi, for example n-2 for xi(n-2)
            out.println(popid+"\t"+nseq+"\t"+L+"\t"+xibegin+"\t"+xiend);
            
            int nx=nseq;
            Random ran=new Random(seed);
            HashMap breaks=new HashMap();
            if(numBreaks<nx-2)  while(breaks.size()<numBreaks)breaks.put(new Integer(ran.nextInt(nx-2)+3), "");//randomly pick numBreaks breaks from 3 to nx
            else for(int i=3;i<=nx;i++)breaks.put(new Integer(i),"");
            numBreaks=breaks.size();
            Integer[] brk=(Integer[])breaks.keySet().toArray(new Integer[numBreaks]);
            Arrays.sort(brk);
            for(int i=0;i<numBreaks;i++)out.print(brk[i]+"\t");
            out.println();//out.println(in.readLine()+"\t"+L);
            if(ii%1==0)System.out.println(ii);
            out.println();
            out.println();
            out.println();
            line=in.readLine();
            out.println(line);
            t=new StringTokenizer(line);
            double[] c=new double[nseq];//count xi(i)
            double total=0;
            for(int i=1;i<=nseq-1;i++) {
                c[i]=Double.parseDouble(t.nextToken());
                total+=c[i];
            }
            c[0]=L-total;
            //create testing set
            double[] ctesting=new double[nseq];
            double Ltesting=(L*(1-ptraining));
            double L0=L;
            
            for(int i=0;i<Ltesting;i++){
                double tmp=(ran.nextDouble()*L);
                for(int j=0;j<nseq;j++){
                    if(tmp>c[j])tmp-=c[j];
                    else{
                        c[j]--;
                        ctesting[j]++;
                        L--;
                        break;
                    }
                }
            }
            SFS_log_likelihood_problem_no_dim_penalty_unfold problem = new SFS_log_likelihood_problem_no_dim_penalty_unfold(nx);
            boolean[] obs=new boolean[nx];
            Arrays.fill(obs,false);
            for(int i=xibegin;i<=xiend;i++) obs[i]=true;
            obs[0]=true;//not necessory
            problem.setObsXi(obs);
            
            Optimizer optimizer = new DE(problem);
            double[] parameters = { 61.9887, 0.6254, 0.4677 };
            boolean keepgoing=true;
            boolean[] splitbefore=new boolean[nx+1];//split before theta i
            
            
            
            splitbefore[2]=true;// dim=1, all thetas are equal, no split
            int dim=1;
            
            int[][] group=new int[dim][];
            int current=0;
            int count=1;
            for(int i=3;i<=nx;i++){
                if(!splitbefore[i])count++;
                else {
                    group[current]=new int[count];
                    count=1;
                    current++;
                }
            }
            group[current]=new int[count];
            current=0;
            group[0][0]=2;
            count=1;
            for(int i=3;i<=nx;i++){
                if(!splitbefore[i]){
                    group[current][count]=i;
                    count++;
                }
                else {
                    current++;
                    group[current][0]=i;
                    count=1;
                }
            }
            
            int numIterations = dimFactor * dim;
            problem.setThetaGroup(group);
            problem.setData(c);
            if(pautocorr<1)problem.setAutoCorr(pautocorr);
            Globals.random = new swarmops.random.MersenneTwister(seed);
            problem.maxIterations = numIterations;
            SFS_log_likelihood_problem_no_dim_penalty_unfold problem2 = new SFS_log_likelihood_problem_no_dim_penalty_unfold(nx);
            problem2.setObsXi(obs);
            problem2.setData(ctesting);
                           
            double bestfit=Double.MAX_VALUE;
            double[] bestest=new double[dim];
            double bestfit2=Double.MAX_VALUE;
            
            int nfeasible=0;
            for(int i=0;i<numRuns;i++){
                Result result=optimizer.optimize(parameters);
                totalRuns++;
                if(result.feasible){
                    //get fitness for testing set
                    problem2.setThetaGroup(group);
                    double fitness=problem2.fitness(result.parameters);
                    if(result.fitness<bestfit){
                        bestfit=result.fitness;
                        bestfit2=fitness;
                        bestest=new double[dim];
                        System.arraycopy(result.parameters, 0, bestest, 0, dim);
                        if(verbose){
                            System.out.println(i+" "+result.iterations+" "+result.fitness+" "+fitness+" "+result.feasible+" "+result.parameters.length);
                            for(int j=0;j<dim;j++)System.out.print(bestest[j]+" ");
                            System.out.println();
                        }
                    }
                    nfeasible++;
                }
                else i--;
            }
            boolean[] currentsplit=(boolean[])splitbefore.clone();
            boolean[] bestsplit=(boolean[])splitbefore.clone();
            if(fitnesstrail){
                out.println("fitness trail:");
                out.print("dim:\t"+dim+"\tbestfit(testing):\t"+bestfit2+"\tbestfit(training):\t"+bestfit);
                for(int i=2;i<=nx;i++){
                    if(bestsplit[i])out.print("\t"+i);
                    else out.print(","+i);
                }
                for(int i=0;i<dim;i++) out.print("\t"+bestest[i]);
                out.println();
            }
            int theta_begin=3;
            while(keepgoing){
                double previous_bestfit=bestfit;
                dim++;
                numIterations = dimFactor * dim;
                problem.maxIterations = numIterations;
                keepgoing=false;
                currentsplit=(boolean[])bestsplit.clone();
                for(int j=theta_begin;j<=nx;j++){
                    splitbefore=(boolean[])currentsplit.clone();
                    if(!splitbefore[j]&&breaks.containsKey(new Integer(j)))splitbefore[j]=true;//only if not in the current split and caintined in the randomly picked breaks
                    else continue;
                    group=new int[dim][];
                    current=0;
                    count=1;
                    for(int i=3;i<=nx;i++){
                        if(!splitbefore[i])count++;
                        else {
                            group[current]=new int[count];
                            count=1;
                            current++;
                        }
                    }
                    group[current]=new int[count];
                    current=0;
                    group[0][0]=2;
                    count=1;
                    for(int i=3;i<=nx;i++){
                        if(!splitbefore[i]){
                            group[current][count]=i;
                            count++;
                        }
                        else {
                            current++;
                            group[current][0]=i;
                            count=1;
                        }
                    }
                    problem.setThetaGroup(group);
                    nfeasible=0;
                    for(int i=0;i<numRuns;i++){
                        Result result=optimizer.optimize(parameters);
                        totalRuns++;
                        if(result.feasible){
                            //get fitness for testing set
                            problem2.setThetaGroup(group);
                            double fitness=problem2.fitness(result.parameters);
                            if(verbose)System.out.println(dim+" "+j+" "+result.fitness);
                            if(result.fitness<bestfit){
                                boolean theta0=false;
                                for(int k=0;k<dim;k++) if(result.parameters[k]==0) theta0=true;
                                if(!theta0){//avoid any theta=0
                                    //keepgoing=true;
                                    bestfit=result.fitness;
                                    bestfit2=fitness;
                                    bestest=new double[dim];
                                    System.arraycopy(result.parameters, 0, bestest, 0, dim);
                                    bestsplit=(boolean[])splitbefore.clone();
                                    if(verbose){
                                        System.out.println(i+" "+result.iterations+" "+result.fitness+" "+fitness+" "+result.feasible+" "+result.parameters.length);
                                        for(int jj=0;jj<dim;jj++)System.out.print(bestest[jj]+" ");
                                        System.out.println();
                                    }
                                }
                            }
                            nfeasible++;
                        }
                        else i--;
                    }
                }
                if(bestfit+3.841/2<previous_bestfit)keepgoing=true;//Chi square test with 1 df, alpha=0.05
                
                if(fitnesstrail&&keepgoing){
                    
                    out.print("dim:\t"+dim+"\tbestfit(testing):\t"+bestfit2+"\tbestfit(training):\t"+bestfit);
                    for(int i=2;i<=nx;i++){
                        if(bestsplit[i])out.print("\t"+i);
                        else out.print(","+i);
                    }
                    for(int i=0;i<dim;i++) out.print("\t"+bestest[i]);
                    out.println();
                }
            }
            dim=bestest.length;
            double[] solution=new double[dim];
            out.println("final model:\t"+bestfit2+"\t"+bestfit);
            for(int i=0;i<dim;i++)solution[i]=bestest[i]*L0;
            for(int i=2;i<=nx;i++){
                if(bestsplit[i])out.print("\t"+i);
                else out.print(","+i);
            }
            out.println();
            for(int i=0;i<dim;i++) out.print(" "+solution[i]);
            out.println();

        }
        
        in.close();
        out.close();
    }

}
