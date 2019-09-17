package edu.mcw.rgd;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.process.Utils;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

public class Manager {

    private DAO dao;
    private String version;
    private String pipelineName;
    private List<String> speciesProcessed;
    private Map< MappedOrtholog,Integer> indexedGenes;
    private  Map<MappedOrtholog,Integer> indexedOrthologs;

    Logger log = Logger.getLogger("core");

    Logger logSummary = Logger.getLogger("status");

    public static void main(String[] args) throws Exception {

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Manager manager = (Manager) (bf.getBean("manager"));

        List<String> species = manager.getSpeciesProcessed();
        try {
            for(String species2: species ){
                manager.run(species2);
            }

        } catch (Exception e) {
            manager.log.error(e);
            throw e;
        }
    }

    void run(String speciesName) throws Exception {

        long startTime = System.currentTimeMillis();

        int speciesTypeKey1 = SpeciesType.parse(speciesName);
        if( SpeciesType.getTaxonomicId(speciesTypeKey1)==0 ) {
            throw new Exception("ERROR: invalid species1: "+speciesName);
        }

        String species1 = SpeciesType.getCommonName(speciesTypeKey1);
        log.info("START: species = " + species1);
        logSummary.info("Summary for species "+ species1);
        List<String> species = getSpeciesProcessed();
        for(String species2: species) {
            int speciesTypeKey2 = SpeciesType.parse(species2);
            if( SpeciesType.getTaxonomicId(speciesTypeKey2)==0 ) {
                throw new Exception("ERROR: invalid species2: "+species2);
            }

            if(speciesTypeKey1 != speciesTypeKey2) {
                log.info("START: Synteny between species1 = " + species1 + " and species2 = " + species2);
                handle(speciesTypeKey1, speciesTypeKey2);
            }
        }
        log.info("END:  time elapsed: " + Utils.formatElapsedTime(startTime, System.currentTimeMillis()));
        log.info("===");
    }

    void handle(int speciesTypeKey1, int speciesTypeKey2) throws Exception {

        int mapKey1 = dao.getPrimaryAssembly(speciesTypeKey1);
        int mapKey2 = dao.getPrimaryAssembly(speciesTypeKey2);
        List<MappedOrtholog> orthologs = dao.getAllOrthologs(speciesTypeKey1,speciesTypeKey2,mapKey1,mapKey2);

        log.info("Remove gene overlaps for Rat");
        List<MappedOrtholog> scannedList1 = removeRatGeneOverlaps(orthologs);

        log.info("Sort the list by Ortholog chromosome and start pos");
        scannedList1.sort((o1, o2) -> {
            int cmp = o1.getDestChromosome().compareTo(o2.getDestChromosome());
            if (cmp == 0)
                cmp = Long.compare(o1.getDestStartPos(), o2.getDestStartPos());
            return cmp;
        });

        log.info("Remove gene overlaps for orthologs");
        List<MappedOrtholog> result = removeDestGeneOverlaps(scannedList1);

        indexedOrthologs = result.stream().collect(
                Collectors.toMap( x -> x, x -> result.indexOf(x)));

        log.info("Sort the list by Rat chromosome and start pos");
        result.sort((o1, o2) -> {
            int cmp = o1.getSrcChromosome().compareTo(o2.getSrcChromosome());
            if (cmp == 0)
                cmp = Long.compare(o1.getSrcStartPos(), o2.getSrcStartPos());
            return cmp;
        });

         indexedGenes = result.stream().collect(
                Collectors.toMap(x -> x,x -> result.indexOf(x)));

        log.info("Generating the blocks");
        List<SyntenyBlock> blocks = generateBlocks(result);

        //Clean up existing blocks
        dao.deleteSyntenyBlocks(mapKey1,mapKey2);

        for(SyntenyBlock block: blocks) {
            log.info(block.getOrientation() + "," + block.getMappedOrtholog().getSrcChromosome() + ","
                    + block.getMappedOrtholog().getSrcStartPos() + "," + block.getMappedOrtholog().getSrcStopPos() + "," + block.getMappedOrtholog().getDestChromosome() + "," +
                    block.getMappedOrtholog().getDestStartPos() + "," + block.getMappedOrtholog().getDestStopPos());

            dao.insertSynteny(block,mapKey1,mapKey2);
        }

        logSummary.info("Blocks between assemblies "+mapKey1+ " and "+mapKey2+ "=" + blocks.size());


    }

    private boolean merge(MappedOrtholog pair,SyntenyBlock syntenyBlock){

        log.info("Checking merge condition of blocks");
        if(syntenyBlock == null)
            return false;
        int orientation = syntenyBlock.getOrientation();
        int index = indexedOrthologs.get(syntenyBlock.getMappedOrtholog());
        MappedOrtholog ortholog = syntenyBlock.getMappedOrtholog();
        int corientation = 0;
        if(pair.getSrcStrand().equalsIgnoreCase(pair.getDestStrand()))
            corientation = 1;
        else corientation = -1;

        return (pair.getSrcChromosome().equalsIgnoreCase(ortholog.getSrcChromosome()) && pair.getDestChromosome().equalsIgnoreCase(ortholog.getDestChromosome())
            && orientation == corientation &&  indexedOrthologs.get(pair) == index + orientation );

    }

    private SyntenyBlock createBlock(MappedOrtholog pair){

        log.info("Creating the new block");
        SyntenyBlock block = new SyntenyBlock();

        if(pair.getSrcStrand().equalsIgnoreCase(pair.getDestStrand()))
            block.setOrientation(1);
        else block.setOrientation(-1);
        block.setMappedOrtholog(pair);
        return block;
    }

    private SyntenyBlock extendBlock(MappedOrtholog pair, SyntenyBlock block){

        log.info("Extending the block");

        SyntenyBlock extendedblock = new SyntenyBlock();
        extendedblock.setOrientation(block.getOrientation());
        MappedOrtholog extPair;
        MappedOrtholog oldPair = block.getMappedOrtholog();
        extPair = pair;
        extPair.setSrcStartPos(Math.min(oldPair.getSrcStartPos(),pair.getSrcStartPos()));
        extPair.setSrcStopPos(Math.max(oldPair.getSrcStopPos(),pair.getSrcStopPos()));
        extPair.setDestStartPos(Math.min(oldPair.getDestStartPos(),pair.getDestStartPos()));
        extPair.setDestStopPos(Math.max(oldPair.getDestStopPos(),pair.getDestStopPos()));
        extendedblock.index = block.getIndex();
        extendedblock.setMappedOrtholog(extPair);
        return extendedblock;
    }

   private List<SyntenyBlock> generateBlocks(List<MappedOrtholog> orthologs){

       ArrayList<SyntenyBlock> blocks = new ArrayList<>();
       SyntenyBlock currentBlock = null;
       for(MappedOrtholog pair : orthologs){
           if(merge(pair,currentBlock)) {
               blocks.remove(currentBlock);
               currentBlock = extendBlock(pair, currentBlock);
               blocks.add(currentBlock);

           }
           else {
               currentBlock = createBlock(pair);
               blocks.add(currentBlock);
           }


       }

       return (blocks);
   }

    private List<SyntenyBlock> massageBlocks(ArrayList<SyntenyBlock> blocks){


        for(int i = 0; i < blocks.size()-1; i++){
            SyntenyBlock current = blocks.get(i);
            SyntenyBlock next = blocks.get(i+1);

            if(i == 0) {
                MappedOrtholog curortholog = current.getMappedOrtholog();
                curortholog.setSrcStartPos(1);
                current.setMappedOrtholog(curortholog);
                blocks.set(i, current);
            }
            if(!current.getMappedOrtholog().getSrcChromosome().equalsIgnoreCase(next.getMappedOrtholog().getSrcChromosome())) {
                MappedOrtholog nextortholog = next.getMappedOrtholog();
                nextortholog.setSrcStartPos(1);
                next.setMappedOrtholog(nextortholog);
                blocks.set(i+1, next);
            }
            else {
                MappedOrtholog curortholog = current.getMappedOrtholog();
                MappedOrtholog nextortholog = next.getMappedOrtholog();

                long delta = (next.getMappedOrtholog().getSrcStartPos()) - (curortholog.getSrcStopPos());
                long epsilon = 1-delta%2;

                long stopPos = curortholog.getSrcStopPos();
                stopPos += delta/2;
                curortholog.setSrcStopPos(stopPos);
                current.setMappedOrtholog(curortholog);

                long startPos = nextortholog.getSrcStartPos();
                startPos -= (delta/2 - epsilon);
                nextortholog.setSrcStartPos(startPos);
                next.setMappedOrtholog(nextortholog);

                blocks.set(i, current);
                blocks.set(i+1, next);
            }

        }

        blocks.sort((o1, o2) -> {
            int cmp = o1.getMappedOrtholog().getDestChromosome().compareTo(o2.getMappedOrtholog().getDestChromosome());
            if (cmp == 0)
                cmp = Long.compare(o1.getMappedOrtholog().getDestStartPos(), o2.getMappedOrtholog().getDestStartPos());
            return cmp;
        });

        for(int i = 0; i < blocks.size()-1; i++){
            SyntenyBlock current = blocks.get(i);
            SyntenyBlock next = blocks.get(i+1);

            if(i == 0) {
                MappedOrtholog curortholog = current.getMappedOrtholog();
                curortholog.setDestStartPos(1);
                current.setMappedOrtholog(curortholog);
                blocks.set(i, current);
            }
            if(!current.getMappedOrtholog().getDestChromosome().equalsIgnoreCase(next.getMappedOrtholog().getDestChromosome())) {
                MappedOrtholog nextortholog = next.getMappedOrtholog();
                nextortholog.setDestStartPos(1);
                next.setMappedOrtholog(nextortholog);
                blocks.set(i+1, next);
            }
            else {
                MappedOrtholog curortholog = current.getMappedOrtholog();
                MappedOrtholog nextortholog = next.getMappedOrtholog();

                long delta = (next.getMappedOrtholog().getDestStartPos()) - (curortholog.getDestStopPos());
                long epsilon = 1-delta%2;

                long stopPos = curortholog.getDestStopPos();
                stopPos += delta/2;
                curortholog.setDestStopPos(stopPos);
                current.setMappedOrtholog(curortholog);

                long startPos = nextortholog.getDestStartPos();
                startPos -= (delta/2 - epsilon);
                nextortholog.setDestStartPos(startPos);
                next.setMappedOrtholog(nextortholog);

                blocks.set(i, current);
                blocks.set(i+1, next);
            }

        }

        blocks.sort((o1, o2) -> {
            int cmp = o1.getMappedOrtholog().getSrcChromosome().compareTo(o2.getMappedOrtholog().getSrcChromosome());
            if (cmp == 0)
                cmp = Long.compare(o1.getMappedOrtholog().getSrcStartPos(), o2.getMappedOrtholog().getSrcStartPos());
            return cmp;
        });
        return blocks;
    }

    private List<MappedOrtholog> removeRatGeneOverlaps(List<MappedOrtholog> orthologList) throws Exception {

        MappedOrtholog pair = null;
        List<MappedOrtholog> result = new ArrayList<>();
        for (MappedOrtholog row : orthologList) {
            if (pair != null && (pair.getSrcChromosome().equalsIgnoreCase(row.getSrcChromosome()) && pair.getSrcStopPos() > row.getSrcStopPos()))
                continue;
            else {
                result.add(row);
                pair = row;
            }
        }
        return result;
    }

    private List<MappedOrtholog> removeDestGeneOverlaps(List<MappedOrtholog> orthologList) throws Exception {

        MappedOrtholog pair = null;
        List<MappedOrtholog> result = new ArrayList<>();
        for (MappedOrtholog row : orthologList) {
            if (pair != null && (pair.getDestChromosome().equalsIgnoreCase(row.getDestChromosome()) && pair.getDestStopPos() > row.getDestStopPos()))
                continue;
            else {
                result.add(row);
                pair = row;
            }
        }
        return result;
    }

    public DAO getDao() {
        return dao;
    }

    public void setDao(DAO dao) {
        this.dao = dao;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public void setPipelineName(String pipelineName) {
        this.pipelineName = pipelineName;
    }

    public String getPipelineName() {
        return pipelineName;
    }

    public List<String> getSpeciesProcessed() {
        return speciesProcessed;
    }

    public void setSpeciesProcessed(List<String> speciesProcessed) {
        this.speciesProcessed = speciesProcessed;
    }

    class SyntenyBlock implements Comparable{


        private int orientation;
        private int index;
        private MappedOrtholog mappedOrtholog;

        public int getOrientation() {
            return orientation;
        }

        public void setOrientation(int orientation) {
            this.orientation = orientation;
        }

        public int getIndex() {
            return index;
        }

        public void setIndex(int index) {
            this.index = index;
        }

        public MappedOrtholog getMappedOrtholog() {
            return mappedOrtholog;
        }

        public void setMappedOrtholog(MappedOrtholog mappedOrtholog) {
            this.mappedOrtholog = mappedOrtholog;
        }

        @Override
        public int compareTo(Object obj1) {
            SyntenyBlock o1 = (SyntenyBlock)obj1;
            int cmp = (o1.getMappedOrtholog().getSrcChromosome().compareTo(mappedOrtholog.getSrcChromosome()));
            if(cmp == 0)
                cmp = Long.compare(o1.getMappedOrtholog().getSrcStartPos(),mappedOrtholog.getSrcStartPos());
            if(cmp == 0)
                cmp = Long.compare(o1.getMappedOrtholog().getSrcStopPos(),mappedOrtholog.getSrcStopPos());
            if(cmp == 0)
                cmp = (o1.getMappedOrtholog().getDestChromosome().compareTo(mappedOrtholog.getDestChromosome()));
            if(cmp == 0)
                cmp = Long.compare(o1.getMappedOrtholog().getDestStartPos(),mappedOrtholog.getDestStartPos());
            if(cmp == 0)
                cmp = Long.compare(o1.getMappedOrtholog().getDestStopPos(),mappedOrtholog.getDestStopPos());
            if(cmp == 0)
                cmp = Long.compare(o1.getOrientation(),orientation);

            return cmp;
        }
    }
}

