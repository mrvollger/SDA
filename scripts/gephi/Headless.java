import java.util.*;  
import java.awt.Color; 
import java.io.File;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import org.gephi.io.exporter.preview.PNGExporter;
import org.gephi.io.exporter.preview.PDFExporter;
import org.gephi.appearance.plugin.PartitionElementColorTransformer;
import org.gephi.appearance.api.AppearanceController;
import org.gephi.appearance.api.AppearanceModel;
import org.gephi.appearance.api.Ranking;
import org.gephi.appearance.api.*;
import org.gephi.filters.api.FilterController;
import org.gephi.filters.api.Query;
import org.gephi.filters.api.Range;
import org.gephi.filters.plugin.graph.DegreeRangeBuilder.DegreeRangeFilter;
import org.gephi.filters.plugin.graph.EgoBuilder.EgoFilter;
import org.gephi.filters.plugin.operator.INTERSECTIONBuilder.IntersectionOperator;
import org.gephi.filters.plugin.partition.PartitionBuilder.NodePartitionFilter;
import org.gephi.graph.api.DirectedGraph;
import org.gephi.graph.api.UndirectedGraph;
import org.gephi.graph.api.Graph;
import org.gephi.graph.api.GraphController;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.GraphView;
import org.gephi.io.importer.api.Container;
import org.gephi.io.importer.api.ImportController;
import org.gephi.io.processor.plugin.DefaultProcessor;
import org.gephi.project.api.ProjectController;
import org.gephi.project.api.Workspace;
import org.gephi.preview.api.PreviewModel;
import org.gephi.preview.api.PreviewController;
import org.gephi.statistics.plugin.GraphDistance;
import org.gephi.appearance.plugin.UniqueLabelColorTransformer;
import org.gephi.io.exporter.api.ExportController;
import org.gephi.preview.api.PreviewProperty;
import org.gephi.preview.types.EdgeColor;
import org.gephi.graph.api.Node;
import org.gephi.graph.api.Edge;
import org.gephi.appearance.spi.RankingTransformer;
import org.gephi.graph.api.Column;
import org.gephi.appearance.api.Function;
import org.gephi.appearance.api.Partition;
import org.gephi.appearance.api.PartitionFunction;
import org.gephi.appearance.plugin.PartitionElementColorTransformer;
import org.gephi.appearance.plugin.palette.Palette;
import org.gephi.appearance.plugin.palette.PaletteManager;
import org.gephi.appearance.plugin.RankingElementColorTransformer;
import org.gephi.appearance.plugin.RankingLabelSizeTransformer;
import org.gephi.appearance.plugin.RankingNodeSizeTransformer;
import org.gephi.preview.api.PreviewController;
import org.gephi.preview.api.PreviewModel;
import org.gephi.preview.api.PreviewProperty;
import org.gephi.layout.plugin.fruchterman.FruchtermanReingold;
import org.gephi.layout.plugin.forceAtlas2.ForceAtlas2;
import org.gephi.layout.plugin.force.yifanHu.YifanHuLayout;
import org.gephi.layout.plugin.force.*;
import org.openide.util.Lookup;
import org.gephi.io.exporter.plugin.ExporterGML;



public class Headless
{

    public static void main(String[] args) {

		// read in file name if provided
		String myfile = "mi.cuts.gml";	
		String outfile = "mi.cuts.gml";	// prefix for output file
		if( args.length >= 2 ){
			myfile = args[0];
			outfile = args[1];
		}	
		System.out.println("Input file is: " + myfile);
	
		//Init a project - and therefore a workspace
        ProjectController pc = Lookup.getDefault().lookup(ProjectController.class);
        pc.newProject();
        Workspace workspace = pc.getCurrentWorkspace();

        //Get controllers and models
        ImportController importController = Lookup.getDefault().lookup(ImportController.class);
        GraphModel graphModel = Lookup.getDefault().lookup(GraphController.class).getGraphModel();
        PreviewModel model = Lookup.getDefault().lookup(PreviewController.class).getModel();
        FilterController filterController = Lookup.getDefault().lookup(FilterController.class);
        AppearanceController ac = Lookup.getDefault().lookup(AppearanceController.class);
        AppearanceModel am = ac.getModel();
        Graph graph = graphModel.getGraph();

        //Import file
        Container container;
        try {
            File file = new File(myfile);
            container = importController.importFile(file);
        } catch (Exception ex) {
            ex.printStackTrace();
            return;
        }
        //Append imported data to GraphAPI
        importController.process(container, new DefaultProcessor(), workspace);
        //See if graph is well imported
        System.out.println("Nodes: " + graph.getNodeCount());
        System.out.println("Edges: " + graph.getEdgeCount());
		 

		// Temperaroly filer out edges that are repulsion edges
		// repulsions edges are marked by "color 1" as opposed to "color 0"		
		ArrayList<Edge> negativeEdges = new ArrayList<Edge>();  
		for (Edge e : graph.getEdges().toArray() ) {
			int color = (int) (long) e.getAttribute("color");
			if(color == 1){
				negativeEdges.add(e);
				graph.removeEdge(e);
			}
		} 
		
		
		// run layout algos 
		int its = 100;	
		// run FA at low gravity to get well sperated clusters
        ForceAtlas2 layout = new ForceAtlas2(null);
		layout.setGraphModel(graphModel);
        layout.resetPropertiesValues();
        layout.initAlgo();
        System.out.println("Running Layout Algorithum: Force Atlas2");
        for (int i = 0; i < its && layout.canAlgo(); i++) {
            layout.goAlgo();
        }
        layout.endAlgo();
		// run FR at high gravity to then bring the clusters into a smaller space 
        FruchtermanReingold layout2 = new FruchtermanReingold(null);
        layout2.setGravity(10.01); 
		layout2.setGraphModel(graphModel);
        layout2.resetPropertiesValues();
        layout2.initAlgo();
        System.out.println("Running Layout Algorithum: FruchtermanReingold");
        for (int i = 0; i < its*100 && layout2.canAlgo(); i++) {
            layout2.goAlgo();
        }
		layout2.endAlgo();
		System.out.println("Done Running Layout Algorithum");
 
		
        // add CCID column
        Column ccCol = graphModel.getNodeTable().addColumn("CCID", Integer.class);
        for (Node n : graph.getNodes()) {
            Long d_color = (Long) n.getAttribute("color");
            int i_color = d_color.intValue();
            n.setAttribute(ccCol, i_color);
        }
		int maxPos = 0;
		int minPos = 99999999;
		for (Node n : graph.getNodes()) {
				int pos =  Integer.parseInt((String)n.getAttribute("pos"));

				if (pos > maxPos) {
						maxPos = pos;
				}
				if (pos < minPos) {
						minPos = pos;
				}
		}
				
				
        //List node columns
        for (Column col : graphModel.getNodeTable()) {
            System.out.println(col);
        }

        // colors can be applied to partitions which always come from attribute columns
        Column ccIDColumn = graphModel.getNodeTable().getColumn("CCID");
        //System.out.println(colorColumn);
        Function func = am.getNodeFunction(graph, ccIDColumn, PartitionElementColorTransformer.class);
         if(func != null){
            Partition partition = ((PartitionFunction) func).getPartition();
            Palette palette = PaletteManager.getInstance().generatePalette(partition.size());
            partition.setColors(palette.getColors());
            ac.transform(func);
						for (Node n : graph.getNodes()) {
								int pos =  Integer.parseInt((String)n.getAttribute("pos"));
								float relPos = 0.25f+ 0.75f*((float) (pos - minPos) )/(maxPos-minPos);
								//						n.setAttribute(posCol, relPos);
								
								n.setAlpha(relPos);
						}
        }else{
            System.out.println("Cannot get partition");
            for(Node n : graph.getNodes()){
								int pos =  Integer.parseInt((String)n.getAttribute("pos"));
								float relPos = 0.25f+ 0.75f*((float) (pos - minPos) )/(maxPos-minPos);
								//						n.setAttribute(posCol, relPos);

                n.setColor( Color.BLUE );
								//System.out.println(relPos);
								n.setAlpha(relPos);
            }
        }

				//        Column posCol = graphModel.getNodeTable().addColumn("POS", Double.class);
								
        
        //Size by Degree
        Function degreeRanking = am.getNodeFunction(graph, AppearanceModel.GraphFunction.NODE_DEGREE, RankingNodeSizeTransformer.class);
        RankingNodeSizeTransformer degreeTransformer = (RankingNodeSizeTransformer) degreeRanking.getTransformer();
        degreeTransformer.setMinSize(25);
        degreeTransformer.setMaxSize(40);
        ac.transform(degreeRanking); 
        
        // change the node labels 
        for(Node n : graph.getNodes()) {
            int i_color = (Integer) n.getAttribute("CCID");
            String s_color = Integer.toString(i_color);
            n.setLabel( s_color  );
        }



        //Preview
        PreviewModel previewModel = Lookup.getDefault().lookup(PreviewController.class).getModel();
        previewModel.getProperties().putValue(PreviewProperty.SHOW_NODE_LABELS, Boolean.TRUE);
        //previewModel.getProperties().putValue(PreviewProperty.NODE_LABEL_PROPORTIONAL_SIZE, Boolean.FALSE);
        model.getProperties().putValue(PreviewProperty.NODE_LABEL_FONT, model.getProperties().getFontValue(PreviewProperty.NODE_LABEL_FONT).deriveFont(8));

        //model.getProperties().putValue(PreviewProperty.EDGE_COLOR, new EdgeColor(Color.GRAY));
        model.getProperties().putValue(PreviewProperty.EDGE_THICKNESS, new Float(0.01f));
        // for some reaosn the next line make it so it takes longer of the pdf to laod in default ubunut 
        // but it does not really change the files ize
        model.getProperties().putValue(PreviewProperty.EDGE_OPACITY, new Float(50) );
        model.getProperties().putValue(PreviewProperty.NODE_PER_NODE_OPACITY, true );				
        
        
		//Export
		ExportController ec = Lookup.getDefault().lookup(ExportController.class);
		try {
			ec.exportFile(new File(outfile + ".pdf"));
			System.out.println("Export Happened");
		} catch (IOException ex) {
			ex.printStackTrace();
			return;
		}
		
		// export gml, these have positions for x,y for the nodes
		// graph re-add negative edges
		for (Edge e : negativeEdges ) {
			graph.addEdge(e);
		}
		ExportController ec2 = Lookup.getDefault().lookup(ExportController.class);
		ExporterGML gml = (ExporterGML) ec2.getExporter(".gml");
		try {
			ec2.exportFile(new File("extraCCplots/cc.positions.gml"), gml);
			System.out.println("Export of gml Happened");
		} catch (IOException ex) {
			ex.printStackTrace();
		}

	
	//	
	// PLOT ONLY NEGATIVE EDGES
	//
	
		// remove positve edges
		for (Edge e : graph.getEdges().toArray() ) {
			graph.removeEdge(e);
		}
		// no edges plot
		ExportController ecnone = Lookup.getDefault().lookup(ExportController.class);
		try {
			ecnone.exportFile(new File( "extraCCplots/" + outfile + ".noEdges.pdf"));
			System.out.println("Export Happened of just nodes happened");
		} catch (IOException ex) {
			ex.printStackTrace();
			return;
		}   	
		// graph re-add negative edges
		for (Edge e : negativeEdges ) {
			graph.addEdge(e);
		}
		// negative edges
		ExportController ec3 = Lookup.getDefault().lookup(ExportController.class);
		try {
			ec3.exportFile(new File( "extraCCplots/" + outfile + ".negative.pdf"));
			System.out.println("Export Happened of negative edges happened");
		} catch (IOException ex) {
			ex.printStackTrace();
			return;
		}   

		//
		// plot negative individually for eahc gorup
		//
		// make every combination of neggative edges
		Hashtable<String, ArrayList<Edge>> repEdges = new Hashtable<String, ArrayList<Edge>>();
		for (Edge e : negativeEdges ) {
			Node source = e.getSource();
			Node target = e.getTarget();
			Integer one = (Integer) source.getAttribute("CCID");
			Integer two = (Integer) target.getAttribute("CCID");
			String CCID1 = one.toString();
			String CCID2 = two.toString();
			//make CCID1 the smaller
			if( one > two ){
				CCID1 = two.toString();
				CCID2 = one.toString();
			}

			if(repEdges.get(CCID1 + "_" + CCID2 ) == null ){
				repEdges.put(CCID1 + "_" + CCID2, new ArrayList<Edge>());
			}
			if( repEdges.get(CCID1) == null){
				repEdges.put(CCID1, new ArrayList<Edge>());
			}	
			if( repEdges.get(CCID2) == null){
				repEdges.put(CCID2, new ArrayList<Edge>());
			}
			
			repEdges.get(CCID1).add(e);
			repEdges.get(CCID2).add(e);
			repEdges.get(CCID1 + "_" + CCID2).add(e);
			
			if(CCID1 == CCID2){
				System.out.println("Same:" + CCID1);
			}
		}
		
		// export every combhination of negative edges	
		Set<String> keys = repEdges.keySet();
		for(String CCID: keys){
			// remove edges from last iterations
			for (Edge e : graph.getEdges().toArray() ) {
				graph.removeEdge(e);
			}
			// add current edges	
			ArrayList<Edge> curEdges = repEdges.get(CCID);
			Integer numNeg = curEdges.size();
			for (Edge e : curEdges ) {
				graph.addEdge(e);
			}
			
			// export subset of graph
			ExportController ecPart = Lookup.getDefault().lookup(ExportController.class);
			try {
				ecPart.exportFile(new File("extraCCplots/" + outfile + ".negative." + CCID.toString() + ".NumNegE." + numNeg.toString()+".pdf"));
				System.out.println("Export Happened of negative edges happened");
			} catch (IOException ex) {
				ex.printStackTrace();
				return;
			}   
		}




	}
}



