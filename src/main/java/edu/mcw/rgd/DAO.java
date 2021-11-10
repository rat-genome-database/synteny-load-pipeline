package edu.mcw.rgd;

import edu.mcw.rgd.dao.AbstractDAO;
import edu.mcw.rgd.dao.impl.OrthologDAO;
import edu.mcw.rgd.datamodel.*;

import java.util.*;


/**
 * @author hsnalabolu
 * wrapper to handle all DAO code
 */
public class DAO extends AbstractDAO {

    OrthologDAO orthologDAO = new OrthologDAO();


    public String getConnectionInfo() {
        return orthologDAO.getConnectionInfo();
    }

    /**
     * get all active gene orthologs for given pair of species
     * @param speciesTypeKey1 species type key for first species
     * @param speciesTypeKey2 species type key for second species
     * @return List of Ortholog objects
     * @throws Exception when unexpected error in spring framework occurs
     */
    public List<MappedOrtholog> getAllOrthologs(int speciesTypeKey1, int speciesTypeKey2,int mapKey1,int mapKey2) throws Exception {
        return orthologDAO.getAllMappedOrthologs(speciesTypeKey1, speciesTypeKey2,mapKey1,mapKey2);
    }

    public void deleteSyntenyBlocks(int mapKey1, int mapKey2) throws Exception{
        String sql = "DELETE FROM SYNTENY WHERE MAP_KEY1 = "+mapKey1+ " AND MAP_KEY2 = "+mapKey2;
        update(sql);
    }

    public void insertSynteny(Manager.SyntenyBlock block,int mapKey1, int mapKey2) throws Exception{

        String sql = "INSERT INTO SYNTENY(MAP_KEY1,MAP_KEY2,CHROMOSOME1,CHROMOSOME2,START_POS1,STOP_POS1,START_POS2,STOP_POS2,ORIENTATION) VALUES(\n"+
                      mapKey1+","+mapKey2+",'"+block.getMappedOrtholog().getSrcChromosome()+"','"+block.getMappedOrtholog().getDestChromosome()+"',"+block.getMappedOrtholog().getSrcStartPos()+","+
                block.getMappedOrtholog().getSrcStopPos()+","+block.getMappedOrtholog().getDestStartPos()+","+block.getMappedOrtholog().getDestStopPos()+","+
                block.getOrientation()+")";
        update(sql);
    }
}
