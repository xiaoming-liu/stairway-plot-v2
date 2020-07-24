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
public class Stairbuilder {
    public static void main(String[] args)throws Exception{
        if(args.length!=1){
            System.out.println("Usage: java Stairbuilder blueprint_file");
            System.exit(1);
        }
        String blueprint_file=args[0];
        BufferedReader in=new BufferedReader(new FileReader(blueprint_file));
        String batchfile=blueprint_file+".sh";
        if(System.getProperty("os.name").startsWith("Windows"))batchfile=blueprint_file+".bat";
        String popid="";
        int nseq=0;
        double L=0;
        boolean folded=true;
        double[] sfs=new double[0];
        double dimFactor=2000;
        double pct_training=0.67;
        long seed=0;
        int[] rand=new int[0];
        String project_dir="";
        String stairway_plot_dir="";
        int ninput=0;
        String inputstem="";
        int binstart=1;
        int binend=-1;
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
            else if(tmp.equalsIgnoreCase("whether_folded")){
                folded=Boolean.parseBoolean(t.nextToken().trim());
                if(folded&&binend<0)binend=nseq/2;
                else if(binend<0)binend=nseq-1;
            }
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
            
            
            else if(tmp.equalsIgnoreCase("smallest_size_of_SFS_bin_used_for_estimation"))binstart=Integer.parseInt(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("largest_size_of_SFS_bin_used_for_estimation")){
                binend=Integer.parseInt(t.nextToken().trim());
                if(folded&&binend>nseq/2){
                    System.out.println("Error: binend>nseq/2! Assuming this is for unfolded bins ...");
                    if(nseq-binend>binstart)binstart=nseq-binend;
                    binend=nseq/2;
                }
            }
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
            else if(tmp.equalsIgnoreCase("random_seed"))seed=Integer.parseInt(t.nextToken().trim());
        }
        in.close();
        
        Random ran=new Random();
        if(seed!=0)ran=new Random(seed);
        else {
            seed=ran.nextLong();
            ran=new Random(seed);
        }
        //output batch file
        String comment="# ";
        String mv="mv -f";
        if(System.getProperty("os.name").startsWith("Windows")){
            comment="REM ";
            mv="MOVE /y";
        }
        String sep=File.separator;
        if(inputstem.isEmpty())inputstem=popid+"-";
        PrintWriter out=new PrintWriter(new FileWriter(batchfile),true);
        System.out.println("File: "+ batchfile + " created");
        
        if(!confirm_or_create_dir(correctDirectory(project_dir))) System.exit(1);
        if(!confirm_or_create_dir(correctDirectory(project_dir+sep+"input"))) System.exit(1);
        for(int i=1;i<=ninput;i++){
            PrintWriter out2=new PrintWriter(new FileWriter(correctDirectory(project_dir+sep+"input")+inputstem+i));
            if(folded) out2.println(popid+"\t"+nseq+"\t"+L+"\t"+binstart+"\t"+binend);
            else out2.println(popid+"\t"+nseq+"\t"+L+"\t"+binstart+"\t"+binend);
            for(int j=0;j<sfs.length;j++)out2.print("\t"+sfs[j]);
            out2.println();
            out2.close();
        }
        out.println(comment+"Step 1: create .addTheta files. random_seed="+seed);
        String date="date";
        if(System.getProperty("os.name").startsWith("Windows"))date="echo %date% %time%";
        out.println(date);
        String program="Stairway_unfold_training_testing7";
        if(folded) program="Stairway_fold_training_testing7";
        for(int ii=0;ii<rand.length;ii++){
            for(int i=1;i<=ninput;i++) {
                out.println("java -Xmx1g -cp "+correctDirectory(stairway_plot_dir)+System.getProperty("path.separator")+correctDirectory(stairway_plot_dir)+"swarmops.jar"+" "+program+" "+correctDirectory(project_dir+sep+"input")+inputstem+i+" "+rand[ii]+" "+pct_training+" "+ran.nextLong());
                
            }
        }
        out.println(date);
        out.println(comment+"Step 2: determine number of break points");
        for(int ii=0;ii<rand.length;ii++){
            if(!confirm_or_create_dir(correctDirectory(project_dir+sep+"rand"+rand[ii]))) System.exit(1);
            for(int i=1;i<=ninput;i++) {
                out.println(mv+" "+correctDirectory(project_dir+sep+"input")+inputstem+i+"."+rand[ii]+"_"+pct_training+".addTheta"+" "+correctDirectory(project_dir+sep+"rand"+rand[ii]));
            }
        }
        out.println("java -Xmx1g -cp "+correctDirectory(stairway_plot_dir)+" "+"Stairpainter"+" "+blueprint_file);
        if(System.getProperty("os.name").startsWith("Windows"))out.println(blueprint_file+".plot.bat");
        else out.println("bash "+blueprint_file+".plot.sh");
        out.println(date);
        out.close();
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
                System.out.println("Directory: "+ dir + " created");
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
}
