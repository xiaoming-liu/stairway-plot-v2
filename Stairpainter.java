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

public class Stairpainter {
    public static void main(String[] args)throws Exception{
        if(args.length!=1){
            System.out.println("Usage: java Stairpainter blueprint_file");
            System.exit(1);
        }
        String blueprint_file=args[0];
        BufferedReader in=new BufferedReader(new FileReader(blueprint_file));
        String batchfile=blueprint_file+".plot.sh";
        if(System.getProperty("os.name").startsWith("Windows"))batchfile=blueprint_file+".plot.bat";
        String popid="";
        int nseq=0;
        double L=0;
        boolean folded=true;
        double[] sfs=new double[0];
        double dimFactor=2000;
        double pct_training=0.67;
        int[] rand=new int[0];
        String project_dir="";
        String stairway_plot_dir="";
        int ninput=0;
        String inputstem="";
        while(in.ready()){
            String line=in.readLine();
            if(line.indexOf("#")!=-1)line=line.substring(0,line.indexOf("#"));//remove comments
            StringTokenizer t=new StringTokenizer(line,":#");
            int nt=t.countTokens();
            if(nt!=2)continue;
            String tmp=t.nextToken().trim();
            if(tmp.equalsIgnoreCase("popid"))popid=t.nextToken().trim();
            else if(tmp.equalsIgnoreCase("inputstem"))inputstem=t.nextToken().trim();
            else if(tmp.equalsIgnoreCase("nseq"))nseq=Integer.parseInt(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("L"))L=Double.parseDouble(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("whether_folded"))folded=Boolean.parseBoolean(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("SFS")){
                String temp=t.nextToken();
                StringTokenizer t2=new StringTokenizer(temp);
                int nt2=t2.countTokens();
                if(folded){
                    if(nt2!=(nseq/2)){
                        System.out.println("SFS error: shall have nseq/2 numbers");
                        if(nt2==(nseq-1)){
                            System.out.println("Assuming the nseq-1 numbers are unfolded counts");
                            double[] sfs0=new double[nseq];
                            for(int i=1;i<=nseq-1;i++)sfs0[i]=Double.parseDouble(t2.nextToken());
                            sfs=new double[nseq/2];
                            for(int i=1;i<=nseq/2;i++){
                                if(i!=nseq-i){
                                    sfs[i-1]=sfs0[i]+sfs0[nseq-i];
                                }
                                else sfs[i-1]=sfs0[i];
                            }
                        }
                        else System.exit(1);
                    }
                    else{
                        sfs=new double[nseq/2];
                        for(int i=0;i<nseq/2;i++)sfs[i]=Double.parseDouble(t2.nextToken());
                    }    
                }
                else{
                    if(nt2!=(nseq-1)){
                        System.out.println("SFS error: shall have nseq-1 numbers");
                        System.exit(1);
                    }
                    else{
                        sfs=new double[nseq-1];
                        for(int i=0;i<nseq-1;i++)sfs[i]=Double.parseDouble(t2.nextToken());
                    } 
                }
            }
            else if(tmp.equalsIgnoreCase("dimFactor"))dimFactor=Integer.parseInt(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("pct_training"))pct_training=Double.parseDouble(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("nrand")){
                String temp=t.nextToken();
                StringTokenizer t2=new StringTokenizer(temp);
                int nt2=t2.countTokens();
                rand=new int[nt2];
                for(int i=0;i<nt2;i++)rand[i]=Integer.parseInt(t2.nextToken());
                Arrays.sort(rand);
            }
            else if(tmp.equalsIgnoreCase("project_dir"))project_dir=t.nextToken().trim();
            else if(tmp.equalsIgnoreCase("stairway_plot_dir"))stairway_plot_dir=t.nextToken().trim();
            else if(tmp.equalsIgnoreCase("ninput"))ninput=Integer.parseInt(t.nextToken().trim());
        }
        in.close();
        if(inputstem.isEmpty())inputstem=popid+"-";
        int ngroup=rand.length;
        double[] average=new double[ngroup];
        double[] median=new double[ngroup];
        String sep=File.separator;
        for(int i=0;i<ngroup;i++){
            String dir=correctDirectory(project_dir+sep+"rand"+rand[i]);
            if(!confirm_dir(dir)) System.exit(1);
            File folder = new File(dir);
            File[] listOfFiles = folder.listFiles();
            String[] infiles=folder.list();
            ArrayList list=new ArrayList();
            for(int ii=0;ii<infiles.length;ii++){
                if(infiles[ii].endsWith(".addTheta"))list.add(infiles[ii]);
            }
            int n1=list.size();
            
            String[] infiles1=(String[])list.toArray(new String[n1]);
            double[] fitness1=new double[n1];
            double total1=0;
            for(int ii=0;ii<n1;ii++){
                System.out.println(dir+infiles1[ii]);
                in=new BufferedReader(new FileReader(dir+infiles1[ii]));
                String line=in.readLine();
                double fitness=0;
                while(line.indexOf("final model:")==-1)line=in.readLine();
                StringTokenizer t=new StringTokenizer(line,"\t");
                t.nextToken();
                fitness=Double.parseDouble(t.nextToken());
                fitness1[ii]=fitness;
                total1+=fitness;
                in.close();
            }
            Arrays.sort(fitness1);
            average[i]=total1/n1;
            median[i]=fitness1[n1/2];
            
        }
        int bestaverage=0;
        int bestmedian=0;
        double averagefit=Double.MAX_VALUE;
        double medianfit=Double.MAX_VALUE;
        for(int i=0;i<ngroup;i++){
            if(averagefit-average[i]>3.841/2){
                averagefit=average[i];
                bestaverage=i;
            }
            else break;
        }
        for(int i=0;i<ngroup;i++){
            if(medianfit-median[i]>3.841/2){
                medianfit=median[i];
                bestmedian=i;
            }
            else break;
        }
        int bestrand=Math.min(bestaverage,bestmedian);
        int bestbreak=rand[bestrand];
        //output batch file
        String comment="# ";
        String mv="mv -f";
        String cp="cp -f";
        if(System.getProperty("os.name").startsWith("Windows")){
            comment="REM ";
            mv="MOVE /y";
            cp="COPY >nul";
        }
        PrintWriter out=new PrintWriter(new FileWriter(batchfile),true);
        for(int i=0;i<ngroup;i++) {
            out.println(comment+"rand"+rand[i]+": average logLikelihood="+(-average[i])+" median logLikelihood="+(-median[i]));
        }
        out.println(comment+"best number of break points (average logLikelihood):"+rand[bestaverage]);
        out.println(comment+"best number of break points (median logLikelihood):"+rand[bestmedian]);
        System.out.println(comment+"best number of break points (average logLikelihood):"+rand[bestaverage]);
        System.out.println(comment+"best number of break points (median logLikelihood):"+rand[bestmedian]);
        out.println(comment+"Step 1: estimations for the \"final\" result");
        if(!confirm_or_create_dir(correctDirectory(project_dir+sep+"final"))) System.exit(1);
        
        for(int i=1;i<=ninput;i++) {
            out.println(cp+" "+correctDirectory(project_dir+sep+"rand"+rand[bestrand])+inputstem+i+"."+rand[bestrand]+"_"+pct_training+".addTheta"+" "+correctDirectory(project_dir+sep+"final"));
        }
        out.println(comment+"Step 2: create summaries and plots");   
        out.println("java -Xmx4g -cp "+correctDirectory(stairway_plot_dir)+System.getProperty("path.separator")+correctDirectory(stairway_plot_dir)+"gral-core-0.11.jar"+System.getProperty("path.separator")+correctDirectory(stairway_plot_dir)+"VectorGraphics2D-0.9.3.jar"+" "+"Stairway_output_summary_plot2"+" "+blueprint_file);
        for(int ii=0;ii<rand.length;ii++){
            if(!confirm_dir(correctDirectory(project_dir+sep+"rand"+rand[ii]))) System.exit(1);
            out.println("java -Xmx4g -cp "+correctDirectory(stairway_plot_dir)+System.getProperty("path.separator")+correctDirectory(stairway_plot_dir)+"gral-core-0.11.jar"+System.getProperty("path.separator")+correctDirectory(stairway_plot_dir)+"VectorGraphics2D-0.9.3.jar"+" "+"Stairway_output_summary_plot2"+" "+blueprint_file+" rand"+rand[ii]);
        }
        out.close();
    }
    private static String correctDirectory(String dir){
        StringTokenizer t=new StringTokenizer (dir,"/\\\"\'");
        int nt=t.countTokens();
        String sep=File.separator;
        String dir2="";
        for(int i=0;i<nt;i++)dir2+=t.nextToken()+sep;
        if((dir.charAt(0)=='\"'||dir.charAt(0)=='\'')&&(dir.charAt(1)=='/'||dir.charAt(1)=='\\'))dir2=sep+dir2;
        else if(dir.charAt(0)=='/'||dir.charAt(0)=='\\')dir2=sep+dir2;
        return dir2;
    }
    private static boolean confirm_dir(String dir){
        File f = new File(dir);
        if (!f.exists()){
            
                System.out.println("Error: Directory: "+ dir + "does not exists");
                return false;
            
        }
        else if(!f.isDirectory()){
            System.out.println("Error: "+dir + "exists but is not a directory");
            return false;
        }
        return true;
    }
    
    private static boolean confirm_or_create_dir(String dir){
        File f = new File(dir);
        if (!f.exists()){
            boolean success = (new File(dir)).mkdir();
            if (success) {
                System.out.println("Warnning: Directory: "+ dir + " created");
                return true;
            }
            else {
                System.out.println("Error: Directory: "+ dir + "does not exists and cannot be created");
                return false;
            }
        }
        else if(!f.isDirectory()){
            System.out.println("Error: "+dir + "exists but is not a directory");
            return false;
        }
        return true;
    }
}
