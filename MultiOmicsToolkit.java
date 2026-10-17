package org.multiomics;

import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Dependency-free Java command-line utilities for the AI Multi-Omics project.
 * Commands:
 *   validate <metadata.csv> <gene_counts.csv> <transcript_counts.csv> <lipidomics.csv>
 *   summarize <matrix.csv> <output.csv>
 *   isoform-usage <transcript_counts.csv> <transcript_gene_map.csv> <output.csv>
 * CSV matrices must have feature_id as the first column and sample IDs as other columns.
 */
public class MultiOmicsToolkit {
    static List<List<String>> readCsv(Path path) throws IOException {
        List<List<String>> rows = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(path)) {
            String line;
            while ((line = br.readLine()) != null) rows.add(parseCsvLine(line));
        }
        if (rows.isEmpty()) throw new IllegalArgumentException("Empty CSV: " + path);
        return rows;
    }
    static List<String> parseCsvLine(String line) {
        List<String> cells = new ArrayList<>(); StringBuilder cell = new StringBuilder(); boolean quoted = false;
        for (int i=0; i<line.length(); i++) {
            char c=line.charAt(i);
            if (c=='"') {
                if (quoted && i+1<line.length() && line.charAt(i+1)=='"') { cell.append('"'); i++; }
                else quoted=!quoted;
            } else if (c==',' && !quoted) { cells.add(cell.toString()); cell.setLength(0); }
            else cell.append(c);
        }
        cells.add(cell.toString()); return cells;
    }
    static String csv(String s) { return "\"" + s.replace("\"", "\"\"") + "\""; }
    static Map<String,Integer> headerMap(List<String> header) {
        Map<String,Integer> m=new LinkedHashMap<>();
        for(int i=0;i<header.size();i++) { String h=header.get(i).trim(); if(m.put(h,i)!=null) throw new IllegalArgumentException("Duplicate column: "+h); }
        return m;
    }
    static double number(String s) {
        if(s==null || s.trim().isEmpty() || s.equalsIgnoreCase("NA") || s.equalsIgnoreCase("NaN")) return Double.NaN;
        try { return Double.parseDouble(s.trim()); } catch(NumberFormatException e) { return Double.NaN; }
    }
    static void validate(String[] args) throws Exception {
        List<List<String>> md=readCsv(Path.of(args[1])); Map<String,Integer> mh=headerMap(md.get(0));
        if(!mh.containsKey("sample_id") || !mh.containsKey("condition")) throw new IllegalArgumentException("Metadata must contain sample_id,condition columns");
        Set<String> samples=new LinkedHashSet<>();
        for(int r=1;r<md.size();r++) { String s=md.get(r).get(mh.get("sample_id")).trim(); if(s.isEmpty() || !samples.add(s)) throw new IllegalArgumentException("Blank or duplicate sample_id: "+s); }
        for(int i=2;i<args.length;i++) {
            Path f=Path.of(args[i]); List<List<String>> rows=readCsv(f); List<String> h=rows.get(0);
            if(h.isEmpty() || !h.get(0).trim().equals("feature_id")) throw new IllegalArgumentException(f+": first column must be feature_id");
            Set<String> featureIds=new HashSet<>();
            for(String s:samples) if(!h.contains(s)) throw new IllegalArgumentException(f+": missing sample column "+s);
            for(int r=1;r<rows.size();r++) {
                if(rows.get(r).size()!=h.size()) throw new IllegalArgumentException(f+": row "+(r+1)+" has "+rows.get(r).size()+" cells; expected "+h.size());
                String id=rows.get(r).get(0).trim(); if(id.isEmpty() || !featureIds.add(id)) throw new IllegalArgumentException(f+": blank/duplicate feature_id "+id);
                for(int c=1;c<h.size();c++) { double v=number(rows.get(r).get(c)); if(!Double.isNaN(v) && (!Double.isFinite(v) || v<0)) throw new IllegalArgumentException(f+": invalid negative/nonfinite value at row "+(r+1)); }
            }
            System.out.println("VALID "+f+": "+(rows.size()-1)+" features; "+samples.size()+" metadata samples");
        }
        System.out.println("Metadata: "+samples.size()+" samples. Input structure checks passed.");
    }
    static void summarize(String[] args) throws Exception {
        List<List<String>> rows=readCsv(Path.of(args[1])); List<String> header=rows.get(0); headerMap(header);
        if(!header.get(0).trim().equals("feature_id")) throw new IllegalArgumentException("First column must be feature_id");
        try(BufferedWriter w=Files.newBufferedWriter(Path.of(args[2]))) {
            w.write("feature_id,n_valid,missing,mean,sd,cv_percent\n");
            for(int r=1;r<rows.size();r++) {
                List<String> row=rows.get(r); if(row.size()!=header.size()) throw new IllegalArgumentException("Wrong column count at row "+(r+1));
                List<Double> vals=new ArrayList<>(); int missing=0;
                for(int c=1;c<row.size();c++){double v=number(row.get(c)); if(Double.isNaN(v)) missing++; else vals.add(v);}
                double mean=Double.NaN,sd=Double.NaN,cv=Double.NaN;
                if(!vals.isEmpty()){mean=vals.stream().mapToDouble(Double::doubleValue).average().orElse(Double.NaN); double ss=0; for(double v:vals)ss+=(v-mean)*(v-mean); if(vals.size()>1)sd=Math.sqrt(ss/(vals.size()-1)); if(mean!=0 && Double.isFinite(sd))cv=100*sd/Math.abs(mean);}
                w.write(csv(row.get(0))+","+vals.size()+","+missing+","+mean+","+sd+","+cv+"\n");
            }
        }
        System.out.println("Feature summaries written to "+args[2]);
    }
    static void isoformUsage(String[] args) throws Exception {
        List<List<String>> mat=readCsv(Path.of(args[1])); List<String> head=mat.get(0); if(!head.get(0).trim().equals("feature_id")) throw new IllegalArgumentException("Transcript matrix first column must be feature_id");
        Map<String,Integer> sampleCol=headerMap(head); List<List<String>> mapRows=readCsv(Path.of(args[2])); Map<String,Integer> mh=headerMap(mapRows.get(0));
        if(!mh.containsKey("transcript_id")||!mh.containsKey("gene_id")) throw new IllegalArgumentException("Mapping needs transcript_id,gene_id columns");
        Map<String,String> txGene=new HashMap<>(); for(int r=1;r<mapRows.size();r++) txGene.put(mapRows.get(r).get(mh.get("transcript_id")).trim(),mapRows.get(r).get(mh.get("gene_id")).trim());
        Map<String,double[]> values=new LinkedHashMap<>(); Map<String,String> gene=new HashMap<>();
        for(int r=1;r<mat.size();r++){String tx=mat.get(r).get(0).trim(); if(!txGene.containsKey(tx)) continue; double[] v=new double[head.size()-1]; for(int c=1;c<head.size();c++){double n=number(mat.get(r).get(c));v[c-1]=Double.isNaN(n)?0:Math.max(0,n);} values.put(tx,v);gene.put(tx,txGene.get(tx));}
        Map<String,double[]> totals=new HashMap<>(); for(String tx:values.keySet()){String g=gene.get(tx);double[] v=values.get(tx);double[] t=totals.computeIfAbsent(g,k->new double[v.length]);for(int i=0;i<v.length;i++)t[i]+=v[i];}
        try(BufferedWriter w=Files.newBufferedWriter(Path.of(args[3]))){w.write("transcript_id,gene_id,sample_id,abundance,isoform_proportion\n");for(String tx:values.keySet()){double[] v=values.get(tx),t=totals.get(gene.get(tx));for(int i=0;i<v.length;i++){double prop=t[i]>0?v[i]/t[i]:0;w.write(csv(tx)+","+csv(gene.get(tx))+","+csv(head.get(i+1))+","+v[i]+","+prop+"\n");}}}
        System.out.println("Per-sample isoform proportions written to "+args[3]+". This is not a differential isoform-switching test.");
    }
    public static void main(String[] args) {
        try {
            if(args.length==0) { usage(); System.exit(2); }
            switch(args[0]) {
                case "validate": if(args.length!=5) throw new IllegalArgumentException("validate expects metadata and three matrices"); validate(args); break;
                case "summarize": if(args.length!=3) throw new IllegalArgumentException("summarize expects input matrix and output CSV"); summarize(args); break;
                case "isoform-usage": if(args.length!=4) throw new IllegalArgumentException("isoform-usage expects transcript matrix, mapping and output CSV"); isoformUsage(args); break;
                default: usage(); System.exit(2);
            }
        } catch(Exception e) { System.err.println("ERROR: "+e.getMessage()); System.exit(1); }
    }
    static void usage(){System.out.println("MultiOmicsToolkit (Java 17+, no external libraries)\n  validate <metadata.csv> <gene_counts.csv> <transcript_counts.csv> <lipidomics.csv>\n  summarize <matrix.csv> <output.csv>\n  isoform-usage <transcript_counts.csv> <transcript_gene_map.csv> <output.csv>");}
}
