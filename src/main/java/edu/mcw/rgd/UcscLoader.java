package edu.mcw.rgd;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;

public class UcscLoader {

    private DAO dao = new DAO();
    private String version;
    private String downloadPrefix;
    private List<String> netFileList;
    private boolean dropAndReload = true;

    Logger log = Logger.getLogger("status");

    public void run() throws Exception {

        log.info(getVersion());
        log.info("   "+dao.getConnectionInfo());

        List<String> processedAssemblies = new ArrayList<>(getNetFileList());
        Collections.shuffle(processedAssemblies);

        processedAssemblies.parallelStream().forEach( netFileInfo -> {

            String[] cols = netFileInfo.split("\\|");
            int mapKey1 = Integer.parseInt(cols[0].trim());
            int mapKey2 = Integer.parseInt(cols[1].trim());
            String localFileName = "data/"+(cols[2].trim());
            String remoteFileName = getDownloadPrefix()+(cols[3].trim());
            boolean loadScaffolds = false;
            if( cols.length>=5 ) {
                String scaffoldText = cols[4].toLowerCase();
                if( scaffoldText.contains("scaffold") ) {
                    loadScaffolds = true;
                }
            }

            FileDownloader fd = new FileDownloader();
            fd.setExternalFile(remoteFileName);
            fd.setLocalFile(localFileName);

            try {
                String localFile = fd.downloadNew();
                run(mapKey1, mapKey2, localFile, loadScaffolds);
            } catch( Exception e ) {
                throw new RuntimeException(e);
            }
        });
    }

    public void run(int mapKey1, int mapKey2, String fileName, boolean loadScaffolds) throws Exception {

        long time0 = System.currentTimeMillis();

        StringBuffer msgBuf = new StringBuffer();

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        msgBuf.append(fileName+" started at "+sdt.format(new Date(time0))+"\n");


        String mapKeyStr = "   [mapKey1="+mapKey1+", mapKey2="+mapKey2+"]\n";

        int blocks = _getCount("SELECT COUNT(*) FROM synteny_ucsc WHERE map_key1="+mapKey1+" AND map_key2="+mapKey2);
        int gaps = _getCount("SELECT COUNT(*) FROM synteny_ucsc_gaps WHERE map_key1="+mapKey1+" AND map_key2="+mapKey2);

        msgBuf.append("  blocks count:    "+Utils.formatThousands(blocks)+mapKeyStr);
        msgBuf.append("  gaps count:    "+Utils.formatThousands(gaps)+mapKeyStr);


        if( dropAndReload ) {
            int rowCount = dao.update("DELETE FROM synteny_ucsc WHERE map_key1=? AND map_key2=?", mapKey1, mapKey2);
            msgBuf.append("  blocks deleted from database:    " + Utils.formatThousands(rowCount)+mapKeyStr);
            rowCount = dao.update("DELETE FROM synteny_ucsc_gaps WHERE map_key1=? AND map_key2=?", mapKey1, mapKey2);
            msgBuf.append("  gaps deleted from database:    " + Utils.formatThousands(rowCount)+mapKeyStr);

        } else {
            // see if we already loaded the data for this file
            int rowCount = dao.getCount("SELECT count(*) FROM synteny_ucsc WHERE map_key1=? AND map_key2=?", mapKey1, mapKey2);
            if (rowCount > 0) {
                msgBuf.append("  blocks already in database:    " + Utils.formatThousands(rowCount)+"\n");
                msgBuf.append("OK -- time elapsed: " + Utils.formatElapsedTime(time0, System.currentTimeMillis())+"\n");
                msgBuf.append("========");
                log.info(msgBuf.toString());
                return;
            }
        }

        BufferedReader in = Utils.openReader(fileName);

        String sql = "INSERT INTO synteny_ucsc (map_key1,map_key2,chromosome1,chromosome2,start_pos1,start_pos2,stop_pos1,stop_pos2," +
                "strand,chain_type,chain_score,chain_level,pos_len1,pos_len2,original_level) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
        String sql2 = "INSERT INTO synteny_ucsc_gaps (map_key1,map_key2,chromosome1,chromosome2,start_pos1,start_pos2,stop_pos1,stop_pos2," +
                "strand,chain_type,chain_score,chain_level,pos_len1,pos_len2,original_level) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";

        Connection conn = DataSourceFactory.getInstance().getDataSource().getConnection();
        conn.setAutoCommit(false);

        PreparedStatement ps = conn.prepareStatement(sql);
        ps.setInt(1, mapKey1);
        ps.setInt(2, mapKey2);

        PreparedStatement psGaps = conn.prepareStatement(sql2);
        psGaps.setInt(1, mapKey1);
        psGaps.setInt(2, mapKey2);

        int linesProcessed = 0;

        String line;
        while( (line=in.readLine())!=null ) {
            String[] cols = line.trim().split("\t",-1);

            linesProcessed++;
            long chainScore = Long.parseLong(cols[11]);

            String tgtChr = cols[2];
            if( tgtChr.startsWith("chr") ) {
                tgtChr = tgtChr.substring(3);
                if( tgtChr.length()>3 ) {
                    continue; // skip unplaced scaffolds
                }
            } else {
                // scaffold GenBank acc
            }

            String srcChr = cols[6];
            if( srcChr.startsWith("chr") ) {
                srcChr = srcChr.substring(3);
                if( srcChr.length()>3 ) {
                    continue; // skip unplaced scaffolds
                }
            } else {
                // scaffold GenBank acc
            }

            int chainLevel = Integer.parseInt(cols[1]); // original level
            int adjLevel = (chainLevel-1)/2 + 1;

            int tgtStart = Integer.parseInt(cols[3]);
            int tgtEnd = Integer.parseInt(cols[4]);
            int tgtLen = tgtEnd-tgtStart;
            String strand = cols[5];
            int srcStart = Integer.parseInt(cols[7]);
            int srcEnd = Integer.parseInt(cols[8]);
            int srcLen = srcEnd-srcStart;

            String chainType = cols[15];

            if( !chainType.equals("gap") ) {
                ps.setString(3, tgtChr);
                ps.setString(4, srcChr);
                ps.setInt(5, tgtStart);
                ps.setInt(6, srcStart);
                ps.setInt(7, tgtEnd);
                ps.setInt(8, srcEnd);
                ps.setString(9, strand);
                ps.setString(10, chainType);
                ps.setLong(11, chainScore);
                ps.setInt(12, adjLevel);
                ps.setInt(13, tgtLen);
                ps.setInt(14, srcLen);
                ps.setInt(15, chainLevel);
                ps.executeUpdate();
            } else {
                psGaps.setString(3, tgtChr);
                psGaps.setString(4, srcChr);
                psGaps.setInt(5, tgtStart);
                psGaps.setInt(6, srcStart);
                psGaps.setInt(7, tgtEnd);
                psGaps.setInt(8, srcEnd);
                psGaps.setString(9, strand);
                psGaps.setString(10, chainType);
                psGaps.setLong(11, chainScore);
                psGaps.setInt(12, adjLevel);
                psGaps.setInt(13, tgtLen);
                psGaps.setInt(14, srcLen);
                psGaps.setInt(15, chainLevel);
                psGaps.executeUpdate();
            }

            if( linesProcessed%250==0 ) {
                conn.commit();
            }
        }
        conn.commit();
        in.close();
        conn.close();


        blocks = _getCount("SELECT COUNT(*) FROM synteny_ucsc WHERE map_key1="+mapKey1+" AND map_key2="+mapKey2);
        gaps = _getCount("SELECT COUNT(*) FROM synteny_ucsc_gaps WHERE map_key1="+mapKey1+" AND map_key2="+mapKey2);

        msgBuf.append("  lines processed: "+Utils.formatThousands(linesProcessed)+mapKeyStr);
        msgBuf.append("  blocks count:    "+Utils.formatThousands(blocks)+mapKeyStr);
        msgBuf.append("  gaps count:    "+Utils.formatThousands(gaps)+mapKeyStr);
        msgBuf.append("OK -- time elapsed: "+Utils.formatElapsedTime(time0, System.currentTimeMillis())+"\n");
        msgBuf.append("========");
        log.info(msgBuf);
    }

    int _getCount(String sql) throws Exception {

        Connection conn = DataSourceFactory.getInstance().getDataSource().getConnection();
        PreparedStatement ps = conn.prepareStatement(sql);
        ResultSet rs = ps.executeQuery();
        int count = -1;
        if( rs.next() ) {
            count = rs.getInt(1);
        }
        conn.close();
        return count;
    }

    public void loadSynNetFile() throws Exception {

        long time0 = System.currentTimeMillis();

        log.info(getVersion());
        log.info("   "+dao.getConnectionInfo());
        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        log.info("   started at "+sdt.format(new Date(time0)));


        String fname = "/Users/mtutaj/downloads/r/rn7.hg38.syn.net";
        BufferedReader in = Utils.openReader(fname);

        BufferedWriter out372_1 = new BufferedWriter(new FileWriter("rn7_1.bed"));
        BufferedWriter out372_2 = new BufferedWriter(new FileWriter("rn7_2.bed"));
        BufferedWriter out372_3 = new BufferedWriter(new FileWriter("rn7_3.bed"));
        BufferedWriter out372_4 = new BufferedWriter(new FileWriter("rn7_4.bed"));
        BufferedWriter out38_1 = new BufferedWriter(new FileWriter("hg38_1.bed"));
        BufferedWriter out38_2 = new BufferedWriter(new FileWriter("hg38_2.bed"));
        BufferedWriter out38_3 = new BufferedWriter(new FileWriter("hg38_3.bed"));
        BufferedWriter out38_4 = new BufferedWriter(new FileWriter("hg38_4.bed"));

        String sql = "INSERT INTO synteny_ucsc (map_key1,map_key2,chromosome1,chromosome2,start_pos1,start_pos2,stop_pos1,stop_pos2," +
                "orientation,chain_type,chain_score,chain_level) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)";
        Connection conn = DataSourceFactory.getInstance().getDataSource().getConnection();
        PreparedStatement ps = conn.prepareStatement(sql);
        int mapKey1 = 372; // rn7
        int mapKey2 = 38;  // hg38
        ps.setInt(1, mapKey1);
        ps.setInt(2, mapKey2);

        String tgtChr = null;
        String line;
        while( (line=in.readLine())!=null ) {
            if( line.startsWith("net") ) {
                String[] cols = line.split(" ");
                tgtChr = cols[1];
                continue;
            }
            // determine level of data
            int level = 0;
            for( int i=0; i<line.length(); i++ ) {
                if( line.charAt(i)==' ' ) {
                    level++;
                } else {
                    break;
                }
            }
            String[] cols = line.trim().split(" ");
            // skip gaps
            if( cols[0].equals("gap") ) {
                continue;
            }
            int tgtStart = Integer.parseInt(cols[1]);
            int tgtLen = Integer.parseInt(cols[2]);
            String srcChr = cols[3];
            String orientation = cols[4]; // + -
            int srcStart = Integer.parseInt(cols[5]);
            int srcLen = Integer.parseInt(cols[6]);

            int chainScore = 0;
            String chainType = null;
            for( int j=7; j<cols.length; j+=2 ) {
                if( cols[j].equals("score") ) {
                    chainScore = Integer.parseInt(cols[j+1]);
                } else if( cols[j].equals("type") ) {
                    chainType = cols[j+1];
                }
            }

            if( srcChr.startsWith("chr") ) {
                srcChr = srcChr.substring(3);
            }
            if( tgtChr.startsWith("chr") ) {
                tgtChr = tgtChr.substring(3);
            }

            int chainLevel = 1 + level/2;
            BufferedWriter out1, out2;
            switch( chainLevel ) {
                case 1:
                    out2 = out38_1;
                    out1 = out372_1;
                    break;
                case 2:
                    out2 = out38_2;
                    out1 = out372_2;
                    break;
                case 3:
                    out2 = out38_3;
                    out1 = out372_3;
                    break;
                case 4:
                    out2 = out38_4;
                    out1 = out372_4;
                    break;
                default:
                    out1 = out2 = null;
            }

            int srcEnd = srcStart+srcLen-1;
            int tgtEnd = tgtStart+tgtLen-1;
            ps.setString(3, tgtChr);
            ps.setString(4, srcChr);
            ps.setInt(5, tgtStart);
            ps.setInt(6, srcStart);
            ps.setInt(7, tgtEnd);
            ps.setInt(8, srcEnd);
            ps.setString(9, orientation);
            ps.setString(10, chainType);
            ps.setInt(11, chainScore);
            ps.setInt(12, chainLevel);
            ps.executeUpdate();

            if( out1!=null ) {
                out1.write("chr"+tgtChr+"\t"+tgtStart+"\t"+tgtEnd+"\t"+orientation+"\n");
            }
            if( out2!=null ) {
                out2.write("chr"+srcChr+"\t"+srcStart+"\t"+srcEnd+"\t"+orientation+"\n");
            }
        }
        in.close();
        conn.close();

        out38_1.close();
        out38_2.close();
        out38_3.close();
        out38_4.close();
        out372_1.close();
        out372_2.close();
        out372_3.close();
        out372_4.close();

        log.info("OK -- time elapsed: "+Utils.formatElapsedTime(time0, System.currentTimeMillis()));
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public String getDownloadPrefix() {
        return downloadPrefix;
    }

    public void setDownloadPrefix(String downloadPrefix) {
        this.downloadPrefix = downloadPrefix;
    }

    public List<String> getNetFileList() {
        return netFileList;
    }

    public void setNetFileList(List<String> netFileList) {
        this.netFileList = netFileList;
    }

}
