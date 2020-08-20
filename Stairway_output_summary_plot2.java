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


import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.io.plots.*;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.axes.LogarithmicRenderer2D;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.graphics.Insets2D;
import de.erichseifert.gral.graphics.Label;
import de.erichseifert.gral.plots.axes.Axis;
import de.erichseifert.gral.plots.axes.LinearRenderer2D;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.io.*;
import java.io.IOException;
import java.util.*;
public class Stairway_output_summary_plot2 {
    public static void main(String[] args)throws Exception{
        if(args.length!=1&&args.length!=2){
            System.out.println("Usage: java Stairway_output_summary_plot2 blueprint_file <output_dir>");
            System.exit(1);
        }
        String blueprint_file=args[0];
        BufferedReader in=new BufferedReader(new FileReader(blueprint_file));
        String project_dir="";
        String title="";
        String xrange="0,0";
        String yrange="0,0";
        double xspacing=2;
        double yspacing=2;
        int fontsize=14;
        double mu=1.2e-8;
        double year_per_generation=24;
        double missing_rate=0;
        int nseq=0;
        double L=0;
        while(in.ready()){
            String line=in.readLine();
            if(line.indexOf("#")!=-1)line=line.substring(0,line.indexOf("#"));//remove comments
            StringTokenizer t=new StringTokenizer(line,":#");
            int nt=t.countTokens();
            if(nt!=2)continue;
            String tmp=t.nextToken().trim();
            if(tmp.equalsIgnoreCase("project_dir"))project_dir=t.nextToken().trim();
            else if(tmp.equalsIgnoreCase("plot_title"))title=t.nextToken().trim();
            else if(tmp.equalsIgnoreCase("nseq"))nseq=Integer.parseInt(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("L"))L=Double.parseDouble(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("xrange"))xrange=t.nextToken().trim();
            else if(tmp.equalsIgnoreCase("yrange"))yrange=t.nextToken().trim();
            else if(tmp.equalsIgnoreCase("xspacing"))xspacing=Integer.parseInt(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("yspacing"))yspacing=Integer.parseInt(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("fontsize"))fontsize=Integer.parseInt(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("mu"))mu=Double.parseDouble(t.nextToken().trim());
            else if(tmp.equalsIgnoreCase("year_per_generation"))year_per_generation=Double.parseDouble(t.nextToken().trim());
        }
        in.close();
        String outfile =correctDirectory(project_dir)+title+".final.summary";
        String sep=File.separator;
        String dir=correctDirectory(project_dir+sep+"final");
        if(args.length==2){
            dir=correctDirectory(project_dir+sep+args[1]);
            outfile =correctDirectory(project_dir)+title+"."+args[1]+".summary";
            title=title+" "+args[1];
        }
        if(!confirm_dir(dir)) System.exit(1);
        File folder = new File(dir);
        File[] listOfFiles = folder.listFiles();
        String[] infiles=folder.list();
        ArrayList list=new ArrayList();
        for(int ii=0;ii<infiles.length;ii++){
            if(infiles[ii].endsWith(".addTheta"))list.add(infiles[ii]);
        }
        int n1=list.size();
        infiles=(String[])list.toArray(new String[n1]);
        
        double xmin=-1;
        double xmax=-1;
        double ymin=-1;
        double ymax=-1;
        if(!xrange.equals("0,0")){
            StringTokenizer t=new StringTokenizer(xrange,", ");
            xmin=Double.parseDouble(t.nextToken());
            xmax=Double.parseDouble(t.nextToken());
        }
        if(!yrange.equals("0,0")){
            StringTokenizer t=new StringTokenizer(yrange,", ");
            ymin=Double.parseDouble(t.nextToken());
            ymax=Double.parseDouble(t.nextToken());
        }
        
        
        
        
        PrintWriter out=new PrintWriter(new FileWriter(outfile),true);
        
        int nrep=infiles.length;
        System.out.println("nrep="+nrep);
        double[][] allest=new double[nseq-1][nrep];
        
       
        for(int ii=0;ii<nrep;ii++){
            in=new BufferedReader(new FileReader(dir+infiles[ii]));
            System.out.println(ii+"\t"+infiles[ii]);
            String line=in.readLine();
            double fitness=0;
            while(line.indexOf("final model:")==-1){
                
                line=in.readLine();
                
            }
            StringTokenizer t=new StringTokenizer(line,"\t");
            t.nextToken();
            fitness=Double.parseDouble(t.nextToken());
            line=in.readLine();
            t=new StringTokenizer(line);
            int dim=t.countTokens();
            int[][] group=new int[dim][];
            for(int i=0;i<dim;i++){
                String tmp=t.nextToken();
                StringTokenizer t2=new StringTokenizer(tmp,",");
                int nxi=t2.countTokens();
                group[i]=new int[nxi];
                for(int j=0;j<nxi;j++)group[i][j]=Integer.parseInt(t2.nextToken());
                
            }
            
            line=in.readLine();
            t=new StringTokenizer(line);
            for(int i=0;i<dim;i++){
                double theta=Double.parseDouble(t.nextToken());
                for(int j=0;j<group[i].length;j++) {
                    allest[group[i][j]-2][ii]=theta;
                    
                }
            }
            
            in.close();
        }
        
        
        int dim=nseq-1;
        double[][] time2=new double[dim][nrep];//record estimated time for each replication
        double[] time3=new double[dim*nrep];//record all time points
        int count=0;
        for(int j=0;j<nrep;j++) {
            for(int i=0;i<dim;i++)time2[i][j]=allest[i][j]/L/(i+2)/(i+1);
            time3[count]=time2[dim-1][j];
            count++;
            for(int i=dim-2;i>=0;i--) {
                time2[i][j]=time2[i][j]+time2[i+1][j];//cummulative time
                time3[count]=time2[i][j];
                count++;
            }
        }
        
        System.out.println("data finished dim="+dim);
        DataTable median=new DataTable(Double.class,Double.class);
        DataTable P125=new DataTable(Double.class,Double.class);
        DataTable P875=new DataTable(Double.class,Double.class);
        DataTable P025=new DataTable(Double.class,Double.class);
        DataTable P975=new DataTable(Double.class,Double.class);
        DataTable median2=new DataTable(Double.class,Double.class);
        float[][] work=new float[dim*nrep][nrep];//record theta at each time points for each replication
        for(int k=0;k<dim*nrep;k++)Arrays.fill(work[k],-1);
        Arrays.sort(time3);
        for(int j=0;j<nrep;j++){
            //System.out.println(j);
            int kbegin=0;
            for(int i=dim-1;i>=0;i--){
                for(int k=kbegin;k<dim*nrep;k++){
                    if(work[k][j]<0&&time2[i][j]>time3[k])work[k][j]=(float)allest[i][j];
                    else if(work[k][j]<0){
                        kbegin=k;
                        break;
                    }
                }
                
            }
            
        }
        out.println("mutation_per_site\tn_estimation\ttheta_per_site_median\ttheta_per_site_2.5%\ttheta_per_site_97.5%\tyear\tNe_median\tNe_2.5%\tNe_97.5%\tNe_12.5%\tNe_87.5%");
        
        
        for(int k=0;k<=dim*nrep-1;k++){
            Arrays.sort(work[k]);
            int npos=0;
            
            for(int i=0;i<nrep;i++){
                if(work[k][i]>0)npos++;
                
            }
            if(npos<nrep*(1-missing_rate))break;
            double[] tmp=new double[npos];
            
            int pos=0;
            int pos2=0;
            int pos3=0;
            for(int i=0;i<nrep;i++){
                if(work[k][i]>0){
                    tmp[pos]=work[k][i];
                    pos++;
                }
                
            }
            if(k==dim*nrep-1){
                
            }
            else{
                
                if(k==0){
                    median.add((1e-100/mu*year_per_generation)/1000,(tmp[npos/2]/L/mu/4)/1000);
                    median2.add((1e-100/mu)/1000,(tmp[npos/2]/L/mu/4)/1000);
                    
                    P025.add((1e-100/mu*year_per_generation)/1000,(tmp[(int)(npos*0.025)]/L/mu/4)/1000);
                    P975.add((1e-100/mu*year_per_generation)/1000,(tmp[(int)(npos*0.975)]/L/mu/4)/1000);
                    P125.add((1e-100/mu*year_per_generation)/1000,(tmp[(int)(npos*0.125)]/L/mu/4)/1000);
                    P875.add((1e-100/mu*year_per_generation)/1000,(tmp[(int)(npos*0.875)]/L/mu/4)/1000);
                    out.println((1e-100)+"\t"+npos+"\t"+(tmp[npos/2]/L)+"\t"+(tmp[(int)(npos*0.025)]/L)+"\t"+(tmp[(int)(npos*0.975)]/L)+"\t"+1e-100/mu*year_per_generation+"\t"+(tmp[npos/2]/L/mu/4)+"\t"+(tmp[(int)(npos*0.025)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.975)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.125)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.875)]/L/mu/4));
                }
                median.add((time3[k]/mu*year_per_generation)/1000,(tmp[npos/2]/L/mu/4)/1000);
                median.add((time3[k+1]/mu*year_per_generation)/1000,(tmp[npos/2]/L/mu/4)/1000);
                median2.add((time3[k]/mu)/1000,(tmp[npos/2]/L/mu/4)/1000);
                median2.add((time3[k+1]/mu)/1000,(tmp[npos/2]/L/mu/4)/1000);
                P025.add((time3[k]/mu*year_per_generation)/1000,(tmp[(int)(npos*0.025)]/L/mu/4)/1000);
                P025.add((time3[k+1]/mu*year_per_generation)/1000,(tmp[(int)(npos*0.025)]/L/mu/4)/1000);
                P975.add((time3[k]/mu*year_per_generation)/1000,(tmp[(int)(npos*0.975)]/L/mu/4)/1000);
                P975.add((time3[k+1]/mu*year_per_generation)/1000,(tmp[(int)(npos*0.975)]/L/mu/4)/1000);
                P125.add((time3[k]/mu*year_per_generation)/1000,(tmp[(int)(npos*0.125)]/L/mu/4)/1000);
                P125.add((time3[k+1]/mu*year_per_generation)/1000,(tmp[(int)(npos*0.125)]/L/mu/4)/1000);
                P875.add((time3[k]/mu*year_per_generation)/1000,(tmp[(int)(npos*0.875)]/L/mu/4)/1000);
                P875.add((time3[k+1]/mu*year_per_generation)/1000,(tmp[(int)(npos*0.875)]/L/mu/4)/1000);
                out.println((time3[k])+"\t"+npos+"\t"+(tmp[npos/2]/L)+"\t"+(tmp[(int)(npos*0.025)]/L)+"\t"+(tmp[(int)(npos*0.975)]/L)+"\t"+time3[k]/mu*year_per_generation+"\t"+(tmp[npos/2]/L/mu/4)+"\t"+(tmp[(int)(npos*0.025)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.975)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.125)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.875)]/L/mu/4));
                out.println((time3[k+1])+"\t"+npos+"\t"+(tmp[npos/2]/L)+"\t"+(tmp[(int)(npos*0.025)]/L)+"\t"+(tmp[(int)(npos*0.975)]/L)+"\t"+time3[k+1]/mu*year_per_generation+"\t"+(tmp[npos/2]/L/mu/4)+"\t"+(tmp[(int)(npos*0.025)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.975)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.125)]/L/mu/4)+"\t"+(tmp[(int)(npos*0.875)]/L/mu/4));
                
            }
        }
        out.close();
        System.out.println("write finished");
        work=new float[0][0];
        time2=new double[0][0];
        time3=new double[0];
        allest=new double[0][0];
        XYPlot plot = new XYPlot();
        
        
        
        nrep=0;
        plot.add(nrep,median,true);
        plot.add(nrep+1,P125,true);
        plot.add(nrep+2,P875,true);
        plot.add(nrep+3,P025,true);
        plot.add(nrep+4,P975,true);
        plot.add(nrep+5,median2,true);
        double insetsTop = 80.0,
               insetsLeft = 90,
               insetsBottom = 80.0,
               insetsRight = 40.0;
        plot.setInsets(new Insets2D.Double(insetsTop, insetsLeft, insetsBottom, insetsRight));
                
        plot.getTitle().setText(title);
        plot.getTitle().setFont(new Font("SansSerif", Font.BOLD, fontsize+4));
        plot.getTitle().setAlignmentY(-50);
        // Style the plot area
        plot.getPlotArea().setBorderColor(new Color(0.0f, 0.3f, 1.0f));
        plot.getPlotArea().setBorderStroke(new BasicStroke(2f));
        Axis x2 = new Axis();
        plot.setAxis(XYPlot.AXIS_X2, x2);
        if(!xrange.equals("0,0")){
            plot.getAxis(XYPlot.AXIS_X).setRange(xmin, xmax);
            plot.getAxis(XYPlot.AXIS_X2).setRange(xmin/year_per_generation, xmax/year_per_generation);
        }
        else {
            double xmin2=plot.getAxis(XYPlot.AXIS_X).getMin().doubleValue();
            double xmax2=plot.getAxis(XYPlot.AXIS_X).getMax().doubleValue();
            plot.getAxis(XYPlot.AXIS_X2).setMin(xmin2/year_per_generation);
            plot.getAxis(XYPlot.AXIS_X2).setMax(xmax2/year_per_generation);
        }
        if(!yrange.equals("0,0"))plot.getAxis(XYPlot.AXIS_Y).setRange(ymin, ymax);
        
        // Style data series
        
        LineRenderer lines1 = new DefaultLineRenderer2D();
        lines1.setColor(new Color(Color.BLACK.getRed(),Color.BLACK.getGreen(),Color.BLACK.getBlue(),255));
        LineRenderer lines2 = new DefaultLineRenderer2D();
        
        LineRenderer lines3 = new DefaultLineRenderer2D();
        lines3.setColor(new Color(Color.CYAN.getRed(),Color.CYAN.getGreen(),Color.CYAN.getBlue(),85));
        LineRenderer lines4 = new DefaultLineRenderer2D();
        lines4.setColor(new Color(Color.DARK_GRAY.getRed(),Color.DARK_GRAY.getGreen(),Color.DARK_GRAY.getBlue(),125));
        LineRenderer lines5 = new DefaultLineRenderer2D();
        lines5.setColor(new Color(Color.GRAY.getRed(),Color.GRAY.getGreen(),Color.GRAY.getBlue(),85));
        LineRenderer lines6 = new DefaultLineRenderer2D();
        lines6.setColor(new Color(Color.GREEN.getRed(),Color.GREEN.getGreen(),Color.GREEN.getBlue(),85));
        LineRenderer lines7 = new DefaultLineRenderer2D();
        lines7.setColor(new Color(Color.LIGHT_GRAY.getRed(),Color.LIGHT_GRAY.getGreen(),Color.LIGHT_GRAY.getBlue(),85));
        LineRenderer lines8 = new DefaultLineRenderer2D();
        lines8.setColor(new Color(Color.MAGENTA.getRed(),Color.MAGENTA.getGreen(),Color.MAGENTA.getBlue(),255));
        LineRenderer lines9 = new DefaultLineRenderer2D();
        lines9.setColor(new Color(Color.ORANGE.getRed(),Color.ORANGE.getGreen(),Color.ORANGE.getBlue(),255));
        LineRenderer lines10 = new DefaultLineRenderer2D();
        lines10.setColor(new Color(Color.PINK.getRed(),Color.PINK.getGreen(),Color.PINK.getBlue(),85));
        LineRenderer lines11 = new DefaultLineRenderer2D();
        lines11.setColor(new Color(Color.RED.getRed(),Color.RED.getGreen(),Color.RED.getBlue(),255));
        LineRenderer lines12 = new DefaultLineRenderer2D();
        lines12.setColor(new Color( Color.YELLOW.getRed(), Color.YELLOW.getGreen(), Color.YELLOW.getBlue(),85));
        LineRenderer lines112 = new DefaultLineRenderer2D();
        lines112.setColor(new Color(Color.RED.getRed(),Color.RED.getGreen(),Color.RED.getBlue(),0));
        
        plot.setLineRenderers(median, lines11);//red
        plot.setPointRenderers(median, null);
        plot.setLineRenderers(median2, lines112);//transparent
        plot.setPointRenderers(median2, null);
        plot.setLineRenderers(P125, lines4);//dark grey
        plot.setPointRenderers(P125, null);
        plot.setLineRenderers(P875, lines4);//dark grey
        plot.setPointRenderers(P875, null);
        plot.setLineRenderers(P025, lines5);//grey
        plot.setPointRenderers(P025, null);
        plot.setLineRenderers(P975, lines5);//grey
        plot.setPointRenderers(P975, null);
        // Style axes
        plot.setAxisRenderer(XYPlot.AXIS_X, new LogarithmicRenderer2D());
        plot.setAxisRenderer(XYPlot.AXIS_Y, new LogarithmicRenderer2D());
        
        plot.setAxisRenderer(XYPlot.AXIS_X2, new LogarithmicRenderer2D());
        plot.setMapping(median2, XYPlot.AXIS_X2, XYPlot.AXIS_Y);
        Label Xlabel=new Label("Time (1k years ago)");
        Xlabel.setFont(new Font("SansSerif", Font.BOLD, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_X).setLabel(Xlabel);
        Label Xlabel2=new Label("Time (1k generations ago)");
        Xlabel2.setFont(new Font("SansSerif", Font.BOLD, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_X2).setLabel(Xlabel2);
        Label Ylabel=new Label("Ne (1k individual)");
        Ylabel.setFont(new Font("SansSerif", Font.BOLD, fontsize));
        Ylabel.setRotation(90.0);
        Ylabel.setAlignmentX(0.3);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setLabel(Ylabel);
        plot.getAxisRenderer(XYPlot.AXIS_X).setIntersection(-Double.MAX_VALUE);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setIntersection( -Double.MAX_VALUE);
        plot.getAxisRenderer(XYPlot.AXIS_X2).setIntersection(Double.MAX_VALUE);
        plot.getAxisRenderer(XYPlot.AXIS_X).setTickFont(new Font("SansSerif", Font.PLAIN, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTickFont(new Font("SansSerif", Font.PLAIN, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_X2).setTickFont(new Font("SansSerif", Font.PLAIN, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_X).setTicksAutoSpaced(false);
        plot.getAxisRenderer(XYPlot.AXIS_X).setTickSpacing(xspacing);
        plot.getAxisRenderer(XYPlot.AXIS_X2).setTickSpacing(xspacing);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTicksAutoSpaced(false);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTickSpacing(yspacing);
        plot.getAxisRenderer(XYPlot.AXIS_X).setMinorTicksVisible(false);
        plot.getAxisRenderer(XYPlot.AXIS_X2).setMinorTicksVisible(false);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setMinorTicksVisible(false);
        
        if(true){
            File file = new File(outfile+".png");
            try {
                    DrawableWriter writer = DrawableWriterFactory.getInstance().get("image/png");
                    writer.write(plot, new FileOutputStream(file), 1500, 750);
            } catch (IOException e) {
            }
        }
        if(true){
            File file2 = new File(outfile+".pdf");
            try {
                    DrawableWriter writer = DrawableWriterFactory.getInstance().get("application/pdf");
                    writer.write(plot, new FileOutputStream(file2), 1500, 750);
            } catch (IOException e) {
            }
        }
            
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
