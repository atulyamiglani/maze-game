import tester.*;
import java.util.function.*;

import javalib.impworld.*;
import java.awt.Color;
import javalib.worldimages.*;
import java.util.*;

/*
 * MODES FOR THE GAME:
 * - BREATH FIRST SEARCH -- OnKeyEvent - "b"
 * - DEPTH FIRST SEARCH -- OnKeyEvent - "d"
 * - PLAYING MODE -- OnKeyEvent - "p"
 * 
 * MOVING IN PLAYING MODE:
 * to move in the maze(after pressing the "p" key for player mode)- OnKeyEvent -
 *  "up", "down", "left", "right" 
 * 
 * OTHER OPTIONS:
 * - RESET THE BOARD -- OnKeyEvent - "r"
 * - CHANGE THE MODE (to get out the current mode, 
 * so that the user can choose the mode again) -- OnKeyEvent - "c"
 * - TOGGLE TO MAKE THE VISITED PATH VISIBLE OR NOT VISIBLE -- OnKeyEvent - "t"
 * 
 * OTHER FUNCTIONALITY:
 * 
 * change the verticalness/ horizontal-ness of the board in the Mazeworld Constructor 
 *
 */

// "b" "d" and "p" for chooseState -- check state updated
// "c" -- check movesCount and graph
// "r" -- check all fields different
// "t" -- check visitedPathsVisible properly changed
// "up" "down" "left" "right" -- check onMove runs

// Represents the states of the cells
enum CellState {
  notVisited, visited, onPath, currCell;
}

// Represents a Cell in a Rectangular Graph used in a Maze game
class Cell {


  //Create cell examples in different states


  ArrayList<Edge> outEdges; // edges connecting with this cell
  int row; // row in the rectangular graph
  int col; // column in the rectangular graph
  CellState state;


  Cell(int row, int col) {
    this.outEdges = new ArrayList<Edge>();
    this.row = row;
    this.col = col;
    this.state = CellState.notVisited;
  }

  // gives the color of the cell based on its state
  Color cellColor(boolean visitedCellsVisible) {
    if (this.state == CellState.notVisited) {
      return Color.GRAY;
    }
    else if (this.state == CellState.visited) {
      if (visitedCellsVisible) {
        return Color.CYAN;
      }
      else {
        return Color.GRAY;
      }

    }
    else if (this.state == CellState.onPath) {
      return Color.BLUE;
    }
    else if (this.state == CellState.currCell) {
      return Color.RED;
    }
    else {
      throw new RuntimeException("Invalid CellState");
    }
  }


  // sets the state of this cell to visited
  boolean hasVisited() {
    return this.state == CellState.visited;
  }

  // sets the state of the cell to isOnPath 
  boolean isOnPath() {
    return this.state == CellState.onPath;
  }

  // checks if the cell is at the given row and column 
  boolean cellAtPos(int row, int col) {
    return this.row == row && this.col == col;
  }


  // sets the state of the cell to visited 
  void setVisited() {
    this.state = CellState.visited;
  }

  //sets the state of the cell to currCell 
  void setCurrCell() {
    this.state = CellState.currCell;
  }

  //sets the state of the cell to onPath
  void setOnPath() {
    this.state = CellState.onPath;
  }

  // applies the given function to all the edges of this cell 
  void forEachEdge(Consumer<Edge> func) {
    for (Edge edge : this.outEdges) {
      func.accept(edge);
    }
  }

  // moves the cell in the given graph by the given row and column 
  Cell moveCell(Graph graph, int dRow, int dCol) {
    return graph.getCell(this.row + dRow, this.col + dCol);
  }



  // Draws this cell in a Maze game. Specifically, draws a gray square with
  // black outlines where the cell does not have an edge (representing a wall).
  WorldImage drawCell(boolean visitedCellsVisible) {

    WorldImage bg = new OverlayImage(new RectangleImage(8, 8, 
        OutlineMode.SOLID, this.cellColor(visitedCellsVisible)), 
        new RectangleImage(10, 10, OutlineMode.SOLID, Color.BLACK));
    for (Edge edge : this.outEdges) {
      if (edge.compareEdgeToCell(this).equals("left")) {
        bg = new OverlayOffsetAlign(AlignModeX.LEFT, AlignModeY.MIDDLE, 
            new RectangleImage(2, 8, OutlineMode.SOLID, 
                this.cellColor(visitedCellsVisible)), 0, 0, bg);

      }
      if (edge.compareEdgeToCell(this).equals("right")) {
        bg = new OverlayOffsetAlign(AlignModeX.RIGHT, AlignModeY.MIDDLE, 
            new RectangleImage(2, 8, OutlineMode.SOLID,
                this.cellColor(visitedCellsVisible)), 0, 0, bg);

      }
      if (edge.compareEdgeToCell(this).equals("bottom")) {
        bg = new OverlayOffsetAlign(AlignModeX.CENTER, AlignModeY.BOTTOM, 
            new RectangleImage(8, 2, OutlineMode.SOLID, 
                this.cellColor(visitedCellsVisible)), 0, 0, bg);

      }
      if (edge.compareEdgeToCell(this).equals("top")) {
        bg = new OverlayOffsetAlign(AlignModeX.CENTER, AlignModeY.TOP, 
            new RectangleImage(8, 2, OutlineMode.SOLID, 
                this.cellColor(visitedCellsVisible)), 0, 0, bg);

      }
    }

    return bg;

  }

  // Compares this cell's position to the given cell's position.
  // Specifically, it returns a string representing where this cell
  // is in relation to the given cell. So if this cell is below the
  // given cell, it returns "bottom". Returns "invalid" if
  // cell is not next to each other.
  String compareCells(Cell c) {
    if (c.row == this.row - 1) {
      return "bottom";
    }
    else if (c.row == this.row + 1) {
      return "top";
    } 
    else if (c.col == this.col + 1) {
      return "left";
    }
    else if (c.col == this.col - 1) {
      return "right";
    }
    return "invalid";
  }

  // EFFECT: Adds the given edge to the outEdges field
  void addEdge(Edge e) {
    this.outEdges.add(e);
  }

  // EFFECT: Removes the given edge from the outEdges field
  void removeEdge(Edge e) {
    this.outEdges.remove(e);
  }

  // checks if the given cell is connected to this cell  
  boolean isConnected(Cell other) {
    for (Edge edge : this.outEdges) {
      if (other.outEdges.contains(edge)) {
        return true;
      }
    }
    return false;
  }

  // returns the edge that connects this cell to the given cell
  public Edge connectedEdge(Cell other) {
    for (Edge edge: this.outEdges) {
      if (edge.getOtherCell(this) == other)  {
        return edge;
      }
    }
    return null;
  }

  // sets the state of this cell to notVisited 
  public void setNotVisited() {
    this.state = CellState.notVisited;
  }
}

// Represents an edge with a weight connecting two cells. In a Maze game,
// having an edge connecting two cells represents the two cells
// having a connection without a wall between.
class Edge {

  Cell from;
  Cell to;
  int weight;

  Edge(Cell from, Cell to, int weight) {
    this.from = from;
    this.to = to;
    this.weight = weight; 
  }

  // EFFECT: Removes this edge from both of its connecting cells.
  void removeThis() {

    this.from.removeEdge(this);
    this.to.removeEdge(this);

  }

  // gets the other cell(between to and from) from this edge 
  // throws an error if the other cell is not next to edge 
  Cell getOtherCell(Cell cell) {
    if (cell == this.from) {
      return this.to;
    }
    else if (cell == this.to) {
      return this.from;
    }
    else {
      throw new RuntimeException("cell not next to edge");
    }
  }

  // Returns the position of this edge in comparison to the given
  // cell. Returns "invalid" if not connected in an orthogonal direction.
  String compareEdgeToCell(Cell c) {
    if (c == this.from) {
      return this.to.compareCells(c);
    }
    else {
      return this.from.compareCells(c);

    }
  }

  // Returns the representative of the from cell
  // in the given union/find data structure.
  Cell getFromRep(HashMap<Cell, Cell> reps) {
    return new Utils().getRep(reps, this.from);
  }


  // Returns the representative of the to cell
  // in the given union/find data structure.
  Cell getToRep(HashMap<Cell, Cell> reps) {
    return new Utils().getRep(reps, this.to);
  }

}

// Represents a graph of cells and edges. In a maze game, this
// graph is "rectangular" and starts out with all of the
// edges in place (i.e. no walls). Only after graphToTree()
// is called is this graph converted into a tree-like graph
// which represents a maze with walls.
class Graph {

  ArrayList<ArrayList<Cell>> allVertices;
  ArrayList<Edge> edges;


  // initializes this graph into a rectangular graph with
  // the given width and height with random weights given
  // by the random generator. No walls have been generated
  // yet -- you have to call graphToTree to generate the walls.
  Graph(int width, int height, Random r, double verticalness) {
    this.edges = new ArrayList<Edge>();

    this.allVertices = new ArrayList<ArrayList<Cell>>();
    for (int row = 0; row < height ; row++) {
      ArrayList<Cell> rowCells = new ArrayList<Cell>();
      for (int col = 0; col < width; col++) {
        rowCells.add(new Cell(row, col));
      }
      allVertices.add(rowCells);
    }

    for (int row = 0; row < height; row++) {
      for (int col = 0; col < width - 1 ; col++) {
        Cell cell1 = allVertices.get(row).get(col);
        Cell cell2 = allVertices.get(row).get(col + 1);
        Edge e = new Edge(cell1, cell2, 
            (int) (r.nextInt(10000) * verticalness)); // change it to a random weight
        cell1.addEdge(e);
        cell2.addEdge(e);
        this.edges.add(e);
      }
    }

    for (int row = 0; row < height - 1; row++) {
      for (int col = 0; col < width ; col++) {
        Cell cell1 = allVertices.get(row).get(col);
        Cell cell2 = allVertices.get(row + 1).get(col);
        Edge e = new Edge(cell1,cell2, r.nextInt(10000)); // change it to a random weight
        cell1.addEdge(e);
        cell2.addEdge(e);
        this.edges.add(e);
      }
    }

  }

  Graph(int width, int height, Random r) {
    this(width, height, r, 1);
  }

  // Returns a rectangular graph with the given width and height
  // with all the connections in place. Call graphToTree to
  // generate the walls.
  Graph(int width, int height) {
    this(width, height, new Random());
  }

  // EFFECT: Converts this graph into a tree using kruskal's algorithm.
  // NOTE: this is slightly different from the given kruskal's algorithm
  // in the assignment because it modifies this graph instead of returning a
  // new graph. The fundamental algorithm is still the same -- the only
  // difference is that this implementation removes edges in this graph when a cycle
  // would be created instead of adding edges into a new graph when a cycle would not
  // be created. We decided that our implementation is better since it avoids copying and
  // creating a whole new graph. This is awkward because the old graph isn't necessary anymore
  // and because it may result in bugs since both the old graph and the new graph
  // have the same cell and edge identities.
  void graphToTree() {

    HashMap<Cell, Cell> reps = new HashMap<Cell, Cell>();

    for (int row = 0; row < this.allVertices.size(); row++ ) {
      for (int col = 0; col < this.allVertices.get(row).size(); col++) {
        reps.put(this.allVertices.get(row).get(col), this.allVertices.get(row).get(col));
      }
    }

    ArrayList<Edge> edges = new ArrayList<Edge>(this.edges);

    // NOTE: field of field okay here because it has been allowed in the past
    // when we've made functional interfaces (in this case, inside a lambda)
    edges.sort((Edge e1, Edge e2) -> e1.weight - e2.weight);

    for (Edge edge : edges) {
      if (edge.getFromRep(reps) == edge.getToRep(reps)) {
        edge.removeThis();
        this.edges.remove(edge);
      }
      else {
        reps.put(edge.getFromRep(reps), edge.getToRep(reps));  
      }
    }
  }

  // Gets the cell at the given row and column positions.
  Cell getCell(int row, int col) {
    if (0 <= row && row < this.allVertices.size() && 
        0 <= col && col < this.allVertices.get(0).size()) {
      return this.allVertices.get(row).get(col);
    }

    return null;
  }


  //resets all the cells of this graph to not visited 
  public void reset() {
    for (ArrayList<Cell> row: this.allVertices) {
      for (Cell cell: row) {
        cell.setNotVisited();
      }
    }

  }

  // counts the number of cells in this graph that have been visited 
  int countHasVisited() {
    int count = 0;
    for (ArrayList<Cell> row : this.allVertices) {
      for (Cell cell: row) {
        if (cell.hasVisited() || cell.isOnPath()) {
          count++;          
        }
      }
    }
    return count;
  }

  // counts the number of cells in this graph that have onPath state
  int countOnPath() {
    int count = 0;
    for (ArrayList<Cell> row : this.allVertices) {
      for (Cell cell: row) {
        if (cell.isOnPath()) {
          count++;          
        }
      }
    }
    return count;
  }
}

// Represents the Utils class
class Utils {

  // Gets the representative of the cell in the union/find data structure
  Cell getRep(HashMap<Cell, Cell> reps, Cell cell) {
    if (reps.get(cell) == cell) {
      return cell;
    }
    else {
      return this.getRep(reps, reps.get(cell));
    }

  }
}

// represents the state of the game 
enum State {
  breadthSearch, depthSearch, playerMove, winState, chooseState;
}

/*
 * World with visitedPathVisible false
 * World with player movement in progress.
 * World where next player move will win
 * World with dfs and bfs in progress.
 * World with dfs and bfs and next move wins.
 * World in process of animating win
 * World which finished win
 */


// Represents a Maze game world state
class MazeWorld extends World {
  int width;
  int height;
  Graph graph;
  Cell currCell;
  boolean visitedPathVisible;
  State state;
  HashMap<Cell, Edge> cameFromEdge;
  Queue<Cell> breadthWorkList;
  Stack<Cell> depthWorkList;
  int movesCount;


  // initializes this Maze World with the given width, height, and graph.
  MazeWorld(int width, int height, Graph graph) {
    this.width = width;
    this.height = height;
    this.graph = graph;
    this.currCell = graph.getCell(0, 0);
    this.visitedPathVisible = true;
    this.state = State.chooseState;
    this.cameFromEdge = new HashMap<>();
    this.breadthWorkList = new ArrayDeque<Cell>();
    this.breadthWorkList.add(this.graph.getCell(0, 0));
    this.depthWorkList = new Stack<Cell>();
    this.depthWorkList.add(this.graph.getCell(0, 0));
    this.movesCount = 0;
  }

  // initializes this Maze World with the given width, height, and using
  // the given random generator to generate the walls.
  MazeWorld(int width, int height, Random r, double verticalness) {
    this(width, height, new Graph(width, height, r, verticalness));
    this.graph.graphToTree();
  }

  MazeWorld(int width, int height, Random r) {
    this(width, height, r, 1);
  }

  MazeWorld(int width, int height, double verticalness) {
    this(width, height, new Random(), verticalness);
  }

  // initializes this Maze World with the given width and height using
  // a random generator.
  MazeWorld(int width, int height) {
    this(width, height, new Random());
  }

  // Draws this MazeWorld.
  public WorldScene makeScene() {
    WorldScene world = new WorldScene(10 * width, 10 * height + 100);
    for (int row = 0; row < this.height; row++) {
      for (int col = 0; col < this.width; col++) {
        world.placeImageXY(
            this.graph.getCell(row, col).drawCell(visitedPathVisible), 5 + 10 * col, 5 + 10 * row );
      }
    }

    world.placeImageXY(new TextImage("moves: " + 
        this.movesCount, 15, Color.BLACK), (10 * width) / 2, (10 * height) + 25);

    if (this.state == State.winState) {
      world.placeImageXY(new TextImage(
          (int) (this.graph.countOnPath() * 100.0 / this.graph.countHasVisited()) + 
          " % right", 15, Color.BLACK),
          (10 * width) / 2, (10 * height) + 75);
    }

    return world;
  }

  // transforms the world based on the keys pressed 
  public void onKeyEvent(String key) {


    if (this.state == State.chooseState) {
      if (key.equals("b")) {
        this.state = State.breadthSearch;
      }
      else if (key.equals("d")) {
        this.state = State.depthSearch;
      }
      else if (key.equals("p")) {
        this.state = State.playerMove;
        this.currCell.setCurrCell();
      }
    }

    if (key.equals("c")) {
      this.state = State.chooseState;
      this.graph.reset();
      this.currCell = this.graph.getCell(0, 0);
      this.cameFromEdge = new HashMap<>();
      this.depthWorkList = new Stack<>();
      this.depthWorkList.add(this.graph.getCell(0, 0));
      this.breadthWorkList = new ArrayDeque<>();
      this.breadthWorkList.add(this.graph.getCell(0, 0));
      this.movesCount = 0;
    }
    
    if (key.equals("r")) {
      this.graph = new Graph(this.width, this.height);
      this.graph.graphToTree();
      this.currCell = this.graph.getCell(0, 0);
      this.state = State.chooseState;
      this.cameFromEdge = new HashMap<>();
      this.breadthWorkList = new ArrayDeque<Cell>();
      this.depthWorkList = new Stack<Cell>();
      this.breadthWorkList.add(this.graph.getCell(0, 0));
      this.depthWorkList.add(this.graph.getCell(0, 0));
      this.movesCount = 0;
    }
    
    if (key.equals("t")) {
      this.visitedPathVisible = !this.visitedPathVisible;
    }

    if (this.state == State.playerMove) {
      this.move(key);
    }

  }

  // changes the world based on the keys - up, down, left, right, in player mode
  public void move(String key) {
    if (key.equals("up")) {
      this.moveInDirection(this.currCell.moveCell(this.graph, -1, 0));
    }
    else if (key.equals("down")) {
      this.moveInDirection(this.currCell.moveCell(this.graph, +1, 0));
    }
    else if (key.equals("left")) {
      this.moveInDirection(this.currCell.moveCell(this.graph, 0, -1));
    }
    else if (key.equals("right")) {
      this.moveInDirection(this.currCell.moveCell(this.graph, 0, +1));
    }
  }

  // moves the currCell state based on the given next cell
  public void moveInDirection(Cell next) {

    if (next == null) {
      return;
    }

    if (next.isConnected(currCell)) {
      this.currCell.setVisited();
      this.movesCount++;
      if (!this.cameFromEdge.keySet().contains(next)) {
        this.cameFromEdge.put(next, next.connectedEdge(currCell));
      }
      next.setCurrCell();
      this.currCell = next;
    }

    if (next.cellAtPos(this.height - 1, this.width - 1)) {
      this.state = State.winState;
    }

  }

  // updates the world on each tick 
  public void onTick() {
    if (this.state == State.breadthSearch) {
      this.searchBreadth();
    }
    else if (this.state == State.depthSearch) {
      this.searchDepth();
    }
    else if (this.state == State.winState) {
      this.animateWinPath();
    }
  }

  // searches for the path based on the breadth-first search algorithm 
  void searchBreadth() {
    this.searchHelp(this.breadthWorkList, this.breadthWorkList.poll());
  }

  // searches for the path based on the depth-first search algorithm 
  void searchDepth() {
    this.searchHelp(this.depthWorkList, this.depthWorkList.pop());
  }

  // implements the breadth first and depth first search algorithms 
  void searchHelp(Collection<Cell> workList, Cell next) {
    next.setVisited();
    this.movesCount++;
    if (next.cellAtPos(this.height - 1, this.width - 1)) {
      this.state = State.winState;
      this.currCell = this.graph.getCell(this.height - 1, this.width - 1);
    }
    else {
      next.forEachEdge(edge -> {
        Cell otherCell = edge.getOtherCell(next);
        if (!otherCell.hasVisited()) {
          workList.add(otherCell);
          this.cameFromEdge.put(otherCell, edge);
        }
      });
    }
  }

  // animates the win path using back tracking in player move mode 
  void animateWinPath() {
    this.currCell.setOnPath();
    if (!this.currCell.cellAtPos(0, 0)) {
      this.currCell = this.cameFromEdge.get(this.currCell).getOtherCell(this.currCell);
    }
  }

}

// Represents the examples and tests 
class ExamplesMaze {

  Graph graph1;
  Graph graph2;
  Graph graph3;
  Graph graph4;

  Graph kruskal3;
  Graph kruskal4;

  Graph graph3SearchInProg;
  Graph graph3SearchEnd;
  Graph graph3SolutionProg;
  Graph graph3SolutionFin;
  
  Cell cellNotVisited;
  Cell cellVisited;
  Cell cellOnPath;
  Cell cellCurrCell;

  // prefix indicates open side
  Cell lCell;
  Cell rCell;
  Cell tCell;
  Cell bCell;
  Cell brCell;
  Cell tbCell;
  Cell openCell;
  Cell tlCornerCell;

  // prefix indicates cell position compared to midCell
  Cell midCell;
  Cell leftCell;
  Cell topCell;
  Cell bottomCell;
  Cell rightCell;

  Edge edgeLeft;
  Edge edgeRight;
  Edge edgeTop;
  Edge edgeBottom;

  MazeWorld world3;
  MazeWorld world4;

  MazeWorld playerMove;
  MazeWorld playerMoveFin;
  MazeWorld depthInProg;
  MazeWorld breadthInProg;
  MazeWorld depthFin;
  MazeWorld breadthFin;
  MazeWorld animateWin;
  MazeWorld finWin;

  HashMap<Cell, Cell> reps;

  // initializes the examples
  void initData() {
    this.graph1 = new Graph(3, 3, new Random(1));
    this.graph2 = new Graph(2, 2, new Random(2));

    this.graph3 = new Graph(5, 5, new Random(3));
    this.graph4 = new Graph(10, 10, new Random(4));

    this.kruskal3 = new Graph(5, 5, new Random(3));
    this.kruskal3.graphToTree();

    this.kruskal4 = new Graph(10, 10, new Random(4));
    this.kruskal4.graphToTree();

    this.graph3SearchInProg = new Graph(5, 5, new Random(3));
    this.graph3SearchInProg.graphToTree();
    this.graph3SearchInProg.getCell(0, 0).setVisited();
    this.graph3SearchInProg.getCell(0, 1).setVisited();
    this.graph3SearchInProg.getCell(1, 0).setVisited();
    this.graph3SearchInProg.getCell(2, 0).setVisited();

    this.graph3SearchEnd = new Graph(5, 5, new Random(3));
    this.graph3SearchEnd.graphToTree();
    this.graph3SearchEnd.getCell(0, 0).setVisited();
    this.graph3SearchEnd.getCell(0, 1).setVisited();
    this.graph3SearchEnd.getCell(1, 0).setVisited();
    this.graph3SearchEnd.getCell(2, 0).setVisited();
    this.graph3SearchEnd.getCell(3, 0).setVisited();
    this.graph3SearchEnd.getCell(4, 0).setVisited();
    this.graph3SearchEnd.getCell(4, 1).setVisited();
    this.graph3SearchEnd.getCell(3, 1).setVisited();
    this.graph3SearchEnd.getCell(3, 2).setVisited();
    this.graph3SearchEnd.getCell(3, 3).setVisited();
    this.graph3SearchEnd.getCell(4, 3).setVisited();

    this.graph3SolutionProg = new Graph(5, 5, new Random(3));
    this.graph3SolutionProg.graphToTree();
    this.graph3SolutionProg.getCell(0, 0).setVisited();
    this.graph3SolutionProg.getCell(0, 1).setVisited();
    this.graph3SolutionProg.getCell(1, 0).setVisited();
    this.graph3SolutionProg.getCell(2, 0).setVisited();
    this.graph3SolutionProg.getCell(3, 0).setVisited();
    this.graph3SolutionProg.getCell(4, 0).setVisited();
    this.graph3SolutionProg.getCell(4, 1).setVisited();
    this.graph3SolutionProg.getCell(3, 1).setOnPath();
    this.graph3SolutionProg.getCell(3, 2).setOnPath();
    this.graph3SolutionProg.getCell(3, 3).setOnPath();
    this.graph3SolutionProg.getCell(4, 3).setOnPath();
    this.graph3SolutionProg.getCell(4, 4).setOnPath();

    this.graph3SolutionFin = new Graph(5, 5, new Random(3));
    this.graph3SolutionFin.graphToTree();
    this.graph3SolutionFin.getCell(0, 0).setOnPath();
    this.graph3SolutionFin.getCell(0, 1).setVisited();
    this.graph3SolutionFin.getCell(1, 0).setOnPath();
    this.graph3SolutionFin.getCell(2, 0).setOnPath();
    this.graph3SolutionFin.getCell(3, 0).setOnPath();
    this.graph3SolutionFin.getCell(4, 0).setOnPath();
    this.graph3SolutionFin.getCell(4, 1).setOnPath();
    this.graph3SolutionFin.getCell(3, 1).setOnPath();
    this.graph3SolutionFin.getCell(3, 2).setOnPath();
    this.graph3SolutionFin.getCell(3, 3).setOnPath();
    this.graph3SolutionFin.getCell(4, 3).setOnPath();
    this.graph3SolutionFin.getCell(4, 4).setOnPath();

    lCell = this.kruskal4.getCell(3, 2);
    rCell = this.kruskal4.getCell(7, 8);
    tCell = this.kruskal4.getCell(6, 3);
    bCell = this.kruskal4.getCell(1, 1);
    brCell = this.kruskal4.getCell(5, 3);
    tbCell = this.kruskal4.getCell(7, 7);
    openCell = this.kruskal4.getCell(6, 7);
    tlCornerCell = this.kruskal4.getCell(0, 0);

    midCell = openCell;
    leftCell = this.kruskal4.getCell(6, 6);
    topCell = this.kruskal4.getCell(5, 7);
    bottomCell = this.kruskal4.getCell(7, 7);
    rightCell = this.kruskal4.getCell(6, 8);

    edgeLeft = midCell.outEdges.get(0);
    edgeRight = midCell.outEdges.get(1);
    edgeTop = midCell.outEdges.get(2);
    edgeBottom = midCell.outEdges.get(3);

    cellNotVisited = lCell;
    cellNotVisited.state = CellState.notVisited;
    cellVisited = rCell;
    cellVisited.state = CellState.visited;
    
    cellOnPath = tCell;
    cellOnPath.state = CellState.onPath;
    
    cellCurrCell = bCell;
    cellCurrCell.state = CellState.currCell;

    world3 = new MazeWorld(5, 5, kruskal3);
    world4 = new MazeWorld(10, 10, kruskal4);

    reps = new HashMap<Cell, Cell>();
    reps.put(midCell, midCell);
    reps.put(leftCell, rightCell);
    reps.put(topCell, topCell);
    reps.put(rightCell, topCell);
    reps.put(bottomCell, midCell);

    playerMove = new MazeWorld(5, 5, this.graph3SearchInProg);

    playerMove.cameFromEdge = new HashMap<>();
    playerMove.cameFromEdge.put(playerMove.graph.getCell(1, 0), 
        playerMove.graph.getCell(1, 0).outEdges.get(0));
    playerMove.cameFromEdge.put(playerMove.graph.getCell(0, 1), 
        playerMove.graph.getCell(0, 1).outEdges.get(0));
    playerMove.cameFromEdge.put(playerMove.graph.getCell(2, 0), 
        playerMove.graph.getCell(2, 0).outEdges.get(0));

    playerMove.currCell = playerMove.graph.getCell(2, 0);
    playerMove.state = State.playerMove;
    playerMove.currCell.setCurrCell();

    playerMoveFin = new MazeWorld(5, 5, this.graph3SearchEnd);

    playerMoveFin.cameFromEdge = new HashMap<>();
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(1, 0), 
        playerMoveFin.graph.getCell(1, 0).outEdges.get(0));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(0, 1), 
        playerMoveFin.graph.getCell(0, 1).outEdges.get(0));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(2, 0), 
        playerMoveFin.graph.getCell(2, 0).outEdges.get(0));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(3, 0), 
        playerMoveFin.graph.getCell(3, 0).outEdges.get(0));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(4, 0), 
        playerMoveFin.graph.getCell(4, 0).outEdges.get(1));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(4, 1), 
        playerMoveFin.graph.getCell(4, 1).outEdges.get(0));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(3, 1), 
        playerMoveFin.graph.getCell(3, 1).outEdges.get(1));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(3, 2), 
        playerMoveFin.graph.getCell(3, 2).outEdges.get(0));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(3, 3), 
        playerMoveFin.graph.getCell(3, 3).outEdges.get(0));
    playerMoveFin.cameFromEdge.put(playerMoveFin.graph.getCell(4, 3), 
        playerMoveFin.graph.getCell(4, 3).outEdges.get(1));

    playerMoveFin.currCell = playerMove.graph.getCell(4, 3);
    playerMoveFin.state = State.playerMove;
    playerMoveFin.currCell.setCurrCell();



    depthInProg = new MazeWorld(5, 5, this.graph3SearchInProg);

    depthInProg.cameFromEdge = new HashMap<>();
    depthInProg.cameFromEdge.put(depthInProg.graph.getCell(1, 0), 
        depthInProg.graph.getCell(1, 0).outEdges.get(0));
    depthInProg.cameFromEdge.put(depthInProg.graph.getCell(0, 1), 
        depthInProg.graph.getCell(0, 1).outEdges.get(0));
    depthInProg.cameFromEdge.put(depthInProg.graph.getCell(2, 0), 
        depthInProg.graph.getCell(2, 0).outEdges.get(0));

    depthInProg.state = State.depthSearch;
    depthInProg.depthWorkList = new Stack<>();
    depthInProg.depthWorkList.add(depthInProg.graph.getCell(1, 1));
    depthInProg.depthWorkList.add(depthInProg.graph.getCell(0, 2));
    depthInProg.depthWorkList.add(depthInProg.graph.getCell(3, 0));




    breadthInProg = new MazeWorld(5, 5, this.graph3SearchInProg);

    breadthInProg.cameFromEdge = new HashMap<>();
    breadthInProg.cameFromEdge.put(breadthInProg.graph.getCell(1, 0), 
        breadthInProg.graph.getCell(1, 0).outEdges.get(0));
    breadthInProg.cameFromEdge.put(breadthInProg.graph.getCell(0, 1), 
        breadthInProg.graph.getCell(0, 1).outEdges.get(0));
    breadthInProg.cameFromEdge.put(breadthInProg.graph.getCell(2, 0), 
        breadthInProg.graph.getCell(2, 0).outEdges.get(0));

    breadthInProg.state = State.breadthSearch;
    breadthInProg.breadthWorkList = new ArrayDeque<>();
    breadthInProg.breadthWorkList.add(breadthInProg.graph.getCell(1, 1));
    breadthInProg.breadthWorkList.add(breadthInProg.graph.getCell(0, 2));
    breadthInProg.breadthWorkList.add(breadthInProg.graph.getCell(3, 0));

    depthFin = new MazeWorld(5, 5, this.graph3SearchEnd);

    depthFin.cameFromEdge = new HashMap<>();
    depthFin.cameFromEdge.put(depthFin.graph.getCell(1, 0), 
        depthFin.graph.getCell(1, 0).outEdges.get(0));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(0, 1), 
        depthFin.graph.getCell(0, 1).outEdges.get(0));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(2, 0), 
        depthFin.graph.getCell(2, 0).outEdges.get(0));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(3, 0), 
        depthFin.graph.getCell(3, 0).outEdges.get(0));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(4, 0), 
        depthFin.graph.getCell(4, 0).outEdges.get(1));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(4, 1), 
        depthFin.graph.getCell(4, 1).outEdges.get(0));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(3, 1), 
        depthFin.graph.getCell(3, 1).outEdges.get(1));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(3, 2), 
        depthFin.graph.getCell(3, 2).outEdges.get(0));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(3, 3), 
        depthFin.graph.getCell(3, 3).outEdges.get(0));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(4, 3), 
        depthFin.graph.getCell(4, 3).outEdges.get(1));
    depthFin.cameFromEdge.put(depthFin.graph.getCell(4, 4), 
        depthFin.graph.getCell(4, 4).outEdges.get(0));

    depthFin.state = State.depthSearch;
    depthFin.depthWorkList = new Stack<>();

    depthFin.depthWorkList.add(depthFin.graph.getCell(1, 1));
    depthFin.depthWorkList.add(depthFin.graph.getCell(0, 2));
    depthFin.depthWorkList.add(depthFin.graph.getCell(4, 2));
    depthFin.depthWorkList.add(depthFin.graph.getCell(3, 4));
    depthFin.depthWorkList.add(depthFin.graph.getCell(4, 4));

    breadthFin = new MazeWorld(5, 5, this.graph3SearchEnd);

    breadthFin.cameFromEdge = new HashMap<>();
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(1, 0), 
        breadthFin.graph.getCell(1, 0).outEdges.get(0));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(0, 1), 
        breadthFin.graph.getCell(0, 1).outEdges.get(0));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(2, 0), 
        breadthFin.graph.getCell(2, 0).outEdges.get(0));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(3, 0), 
        breadthFin.graph.getCell(3, 0).outEdges.get(0));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(4, 0), 
        breadthFin.graph.getCell(4, 0).outEdges.get(1));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(4, 1), 
        breadthFin.graph.getCell(4, 1).outEdges.get(0));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(3, 1), 
        breadthFin.graph.getCell(3, 1).outEdges.get(1));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(3, 2), 
        breadthFin.graph.getCell(3, 2).outEdges.get(0));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(3, 3), 
        breadthFin.graph.getCell(3, 3).outEdges.get(0));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(4, 3), 
        breadthFin.graph.getCell(4, 3).outEdges.get(1));
    breadthFin.cameFromEdge.put(breadthFin.graph.getCell(4, 4), 
        breadthFin.graph.getCell(4, 4).outEdges.get(0));

    breadthFin.state = State.breadthSearch;
    breadthFin.breadthWorkList = new ArrayDeque<>();

    breadthFin.breadthWorkList.add(depthFin.graph.getCell(4, 4));
    breadthFin.breadthWorkList.add(depthFin.graph.getCell(1, 1));
    breadthFin.breadthWorkList.add(depthFin.graph.getCell(0, 2));
    breadthFin.breadthWorkList.add(depthFin.graph.getCell(4, 2));
    breadthFin.breadthWorkList.add(depthFin.graph.getCell(3, 4));

    animateWin = new MazeWorld(5, 5, this.graph3SolutionProg);
    animateWin.currCell = animateWin.graph.getCell(4, 1);

    animateWin.cameFromEdge = new HashMap<>();
    animateWin.cameFromEdge.put(animateWin.graph.getCell(1, 0), 
        animateWin.graph.getCell(1, 0).outEdges.get(0));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(0, 1), 
        animateWin.graph.getCell(0, 1).outEdges.get(0));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(2, 0), 
        animateWin.graph.getCell(2, 0).outEdges.get(0));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(3, 0), 
        animateWin.graph.getCell(3, 0).outEdges.get(0));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(4, 0), 
        animateWin.graph.getCell(4, 0).outEdges.get(1));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(4, 1), 
        animateWin.graph.getCell(4, 1).outEdges.get(0));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(3, 1), 
        animateWin.graph.getCell(3, 1).outEdges.get(1));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(3, 2), 
        animateWin.graph.getCell(3, 2).outEdges.get(0));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(3, 3), 
        animateWin.graph.getCell(3, 3).outEdges.get(0));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(4, 3), 
        animateWin.graph.getCell(4, 3).outEdges.get(1));
    animateWin.cameFromEdge.put(animateWin.graph.getCell(4, 4), 
        animateWin.graph.getCell(4, 4).outEdges.get(0));

    animateWin.state = State.winState;

    finWin = new MazeWorld(5, 5, this.graph3SolutionFin);
    finWin.currCell = finWin.graph.getCell(0, 0);

    finWin.cameFromEdge = new HashMap<>();
    finWin.cameFromEdge.put(finWin.graph.getCell(1, 0), 
        finWin.graph.getCell(1, 0).outEdges.get(0));
    finWin.cameFromEdge.put(finWin.graph.getCell(0, 1), 
        finWin.graph.getCell(0, 1).outEdges.get(0));
    finWin.cameFromEdge.put(finWin.graph.getCell(2, 0), 
        finWin.graph.getCell(2, 0).outEdges.get(0));
    finWin.cameFromEdge.put(finWin.graph.getCell(3, 0), 
        finWin.graph.getCell(3, 0).outEdges.get(0));
    finWin.cameFromEdge.put(finWin.graph.getCell(4, 0), 
        finWin.graph.getCell(4, 0).outEdges.get(1));
    finWin.cameFromEdge.put(finWin.graph.getCell(4, 1), 
        finWin.graph.getCell(4, 1).outEdges.get(0));
    finWin.cameFromEdge.put(finWin.graph.getCell(3, 1), 
        finWin.graph.getCell(3, 1).outEdges.get(1));
    finWin.cameFromEdge.put(finWin.graph.getCell(3, 2), 
        finWin.graph.getCell(3, 2).outEdges.get(0));
    finWin.cameFromEdge.put(finWin.graph.getCell(3, 3), 
        finWin.graph.getCell(3, 3).outEdges.get(0));
    finWin.cameFromEdge.put(finWin.graph.getCell(4, 3), 
        finWin.graph.getCell(4, 3).outEdges.get(1));
    finWin.cameFromEdge.put(finWin.graph.getCell(4, 4), 
        finWin.graph.getCell(4, 4).outEdges.get(0));

    finWin.state = State.winState;

  }

  // DRAWS THE MAZE WORLD IN A NEW WINDOW.
  void testBigBang(Tester t) {
    this.initData();

    // CHANGE THIS IF YOU WANT TO CHANGE HOW THE MAZE WORLD IS INITIALIZED
    MazeWorld m1 = new MazeWorld(100, 60);


    int worldWidth = m1.width * 50;
    int worldHeight = m1.height * 50 + 100;
    m1.bigBang(worldWidth, worldHeight, 0.01);

  }

  // Tests the Cell Constructor
  void testCellConstructor(Tester t) {
    this.initData();

    Cell cell1 = new Cell(0, 0);
    Cell cell2 = new Cell(10, 10);
    Cell cell3 = new Cell(3, 6);

    t.checkExpect(cell1.row, 0);
    t.checkExpect(cell1.col, 0);
    t.checkExpect(cell1.outEdges, new ArrayList<Edge>());

    t.checkExpect(cell2.row, 10);
    t.checkExpect(cell2.col, 10);
    t.checkExpect(cell2.outEdges, new ArrayList<Edge>());

    t.checkExpect(cell3.row, 3);
    t.checkExpect(cell3.col, 6);
    t.checkExpect(cell3.outEdges, new ArrayList<Edge>());
  }

  // Tests drawCell
  void testDrawCell(Tester t) {
    this.initData();

    t.checkExpect(this.lCell.drawCell(false), 
        new OverlayOffsetAlign(AlignModeX.LEFT, AlignModeY.MIDDLE, 
        new RectangleImage(2, 8, OutlineMode.SOLID, Color.GRAY), 0, 0,
        new OverlayImage(new RectangleImage(8, 8, OutlineMode.SOLID, Color.gray), 
            new RectangleImage(10, 10, OutlineMode.SOLID, Color.BLACK))));

    this.rCell.state = CellState.currCell;
    t.checkExpect(this.rCell.drawCell(false), 
        new OverlayOffsetAlign(AlignModeX.RIGHT, AlignModeY.MIDDLE, 
            new RectangleImage(2, 8, OutlineMode.SOLID, Color.RED), 0, 0, 
            new OverlayImage(new RectangleImage(8, 8, OutlineMode.SOLID, Color.RED), 
                new RectangleImage(10, 10, OutlineMode.SOLID, Color.BLACK))));

    this.tCell.state = CellState.onPath;
    t.checkExpect(this.tCell.drawCell(false), 
        new OverlayOffsetAlign(AlignModeX.CENTER, AlignModeY.TOP, 
            new RectangleImage(8, 2, OutlineMode.SOLID, Color.BLUE), 0, 0, 
            new OverlayImage(new RectangleImage(8, 8, OutlineMode.SOLID, Color.BLUE), 
                new RectangleImage(10, 10, OutlineMode.SOLID, Color.BLACK))));

    this.bCell.state = CellState.visited;
    t.checkExpect(this.bCell.drawCell(true), 
        new OverlayOffsetAlign(AlignModeX.CENTER, AlignModeY.BOTTOM, 
            new RectangleImage(8, 2, OutlineMode.SOLID, Color.CYAN), 0, 0, 
            new OverlayImage(new RectangleImage(8, 8, OutlineMode.SOLID, Color.CYAN), 
                new RectangleImage(10, 10, OutlineMode.SOLID, Color.BLACK))));

    this.brCell.state = CellState.visited;
    t.checkExpect(this.brCell.drawCell(false), 
        new OverlayOffsetAlign(AlignModeX.CENTER, AlignModeY.BOTTOM, 
            new RectangleImage(8, 2, OutlineMode.SOLID, Color.GRAY), 0, 0, 
            new OverlayOffsetAlign(AlignModeX.RIGHT, AlignModeY.MIDDLE,
                new RectangleImage(2, 8, OutlineMode.SOLID, Color.GRAY), 0, 0, 
                new OverlayImage(new RectangleImage(8, 8, OutlineMode.SOLID, Color.gray), 
                    new RectangleImage(10, 10, OutlineMode.SOLID, Color.BLACK)))));

    t.checkExpect(this.openCell.drawCell(true), 
        new OverlayOffsetAlign(AlignModeX.CENTER, AlignModeY.BOTTOM, 
            new RectangleImage(8, 2, OutlineMode.SOLID, Color.GRAY), 0, 0, 
            new OverlayOffsetAlign(AlignModeX.CENTER, AlignModeY.TOP, 
                new RectangleImage(8, 2, OutlineMode.SOLID, Color.GRAY), 0, 0, 
                new OverlayOffsetAlign(AlignModeX.RIGHT, AlignModeY.MIDDLE, 
                    new RectangleImage(2, 8, OutlineMode.SOLID, Color.GRAY), 0, 0, 
                    new OverlayOffsetAlign(AlignModeX.LEFT, AlignModeY.MIDDLE, 
                        new RectangleImage(2, 8, OutlineMode.SOLID, Color.GRAY), 0, 0, 
                        new OverlayImage(new RectangleImage(8, 8, OutlineMode.SOLID, Color.gray), 
                            new RectangleImage(10, 10, OutlineMode.SOLID, Color.BLACK)))))));

  }

  // Tests compareCells
  void testCompareCells(Tester t) {
    this.initData();
    t.checkExpect(this.rightCell.compareCells(this.midCell), "right");
    t.checkExpect(this.leftCell.compareCells(this.midCell), "left");
    t.checkExpect(this.topCell.compareCells(this.midCell), "top");
    t.checkExpect(this.bottomCell.compareCells(this.midCell), "bottom");
    t.checkExpect(this.lCell.compareCells(this.midCell), "invalid");
  }

  // Tests addEdge
  void testAddEdge(Tester t) {
    this.initData();
    this.bCell.addEdge(edgeBottom);
    t.checkExpect(this.bCell.outEdges.contains(edgeBottom), true);

    this.bottomCell.addEdge(edgeTop);
    t.checkExpect(this.bottomCell.outEdges.contains(edgeTop), true);

    this.brCell.addEdge(edgeLeft);
    t.checkExpect(this.brCell.outEdges.contains(edgeLeft), true);

  }

  // Tests removeEdge
  void testRemoveEdge(Tester t) {
    this.initData();
    this.leftCell.removeEdge(edgeBottom);
    t.checkExpect(this.leftCell.outEdges.contains(edgeBottom), false);

    this.bottomCell.removeEdge(edgeTop);
    t.checkExpect(this.bottomCell.outEdges.contains(edgeTop), false);

    this.rightCell.removeEdge(edgeRight);
    t.checkExpect(this.rightCell.outEdges.contains(edgeRight), false);
  }
  
  // Tests the Edge constructor
  
  void testEdgeConstructor(Tester t) {
    
    this.initData();

    Edge edge1 = new Edge(lCell, rCell, 0);
    Edge edge2 = new Edge(rCell, tCell, 5);
    Edge edge3 = new Edge(bCell, brCell, 10);

    t.checkExpect(edge1.from, lCell);
    t.checkExpect(edge1.to, rCell);
    t.checkExpect(edge1.weight, 0);

    t.checkExpect(edge2.from, rCell);
    t.checkExpect(edge2.to, tCell);
    t.checkExpect(edge2.weight, 5);

    t.checkExpect(edge3.from, bCell);
    t.checkExpect(edge3.to, brCell);
    t.checkExpect(edge3.weight, 10);
  }

  // Tests removeThis
  void testRemoveThis(Tester t) {
    this.initData();

    t.checkExpect(this.leftCell.outEdges.contains(edgeLeft), true);
    t.checkExpect(this.midCell.outEdges.contains(edgeLeft), true);
    t.checkExpect(this.topCell.outEdges.contains(edgeTop), true);
    t.checkExpect(this.midCell.outEdges.contains(edgeTop), true);

    this.edgeLeft.removeThis();
    this.edgeTop.removeThis();

    t.checkExpect(this.leftCell.outEdges.contains(edgeLeft), false);
    t.checkExpect(this.midCell.outEdges.contains(edgeLeft), false);
    t.checkExpect(this.topCell.outEdges.contains(edgeTop), false);
    t.checkExpect(this.midCell.outEdges.contains(edgeTop), false);
  }


  // Tests compareEdgeToCell
  void testCompareEdgeToCell(Tester t) {
    this.initData();
    t.checkExpect(this.edgeLeft.compareEdgeToCell(midCell), "left");
    t.checkExpect(this.edgeBottom.compareEdgeToCell(midCell), "bottom");
    t.checkExpect(this.edgeRight.compareEdgeToCell(midCell), "right");
    t.checkExpect(this.edgeTop.compareEdgeToCell(midCell), "top");
    t.checkExpect(this.edgeLeft.compareEdgeToCell(bCell), "invalid");
  }

  // Tests getFromRep
  void testGetFromRep(Tester t) {
    this.initData();
    t.checkExpect(this.edgeLeft.getFromRep(reps), this.topCell);
    t.checkExpect(this.edgeRight.getFromRep(reps), this.midCell);
    t.checkExpect(this.edgeBottom.getFromRep(reps), this.midCell);
    t.checkExpect(this.edgeTop.getFromRep(reps), this.topCell);
  }

  // Tests getToRep
  void testGetToRep(Tester t) {
    this.initData();
    t.checkExpect(this.edgeLeft.getToRep(reps), this.midCell);
    t.checkExpect(this.edgeRight.getToRep(reps), this.topCell);
    t.checkExpect(this.edgeBottom.getToRep(reps), this.midCell);
    t.checkExpect(this.edgeTop.getToRep(reps), this.midCell);
  }

  // Tests the Graph constructors
  void testGraphConstructors(Tester t) {
    this.initData();

    Cell tlCell = this.graph2.allVertices.get(0).get(0);
    Cell trCell = this.graph2.allVertices.get(0).get(1);
    Cell blCell = this.graph2.allVertices.get(1).get(0);
    Cell brCell = this.graph2.allVertices.get(1).get(1);

    Edge tEdge = tlCell.outEdges.get(0);
    Edge bEdge = blCell.outEdges.get(0);
    Edge lEdge = tlCell.outEdges.get(1);
    Edge rEdge = trCell.outEdges.get(1);

    t.checkExpect(graph2.edges, new ArrayList<Edge>(Arrays.asList(tEdge, bEdge, lEdge, rEdge)));
    t.checkExpect(graph2.allVertices, new ArrayList<ArrayList<Cell>>(
        Arrays.asList(new ArrayList<Cell>(Arrays.asList(tlCell, trCell)),
            new ArrayList<Cell>(Arrays.asList(blCell, brCell)))));


    Graph randGraph = new Graph(20, 15);

    t.checkExpect(randGraph.edges.size(), 565);
    t.checkExpect(randGraph.allVertices.size(), 15);
    t.checkExpect(randGraph.allVertices.get(0).size(), 20);

    Graph verticalGraph = new Graph(5, 7, new Random(1), 100000.0);
    verticalGraph.graphToTree();
    t.checkExpect(verticalGraph.getCell(0, 0).isConnected(verticalGraph.getCell(1, 0)), true);

    Graph horizontalGraph = new Graph(4, 10, new Random(1), 0.00001);
    horizontalGraph.graphToTree();
    t.checkExpect(horizontalGraph.getCell(0, 0).isConnected(horizontalGraph.getCell(0, 1)), true);
  }

  // Tests graphToTree
  void testGraphToTree(Tester t) {
    this.initData();

    Graph graph1Tree = new Graph(3, 3, new Random(1));

    t.checkExpect(graph1, graph1Tree);

    Edge removedEdge1 = graph1Tree.getCell(0, 0).outEdges.get(0);
    Edge removedEdge2 = graph1Tree.getCell(0, 2).outEdges.get(1);
    Edge removedEdge3 = graph1Tree.getCell(1, 1).outEdges.get(3);
    Edge removedEdge4 = graph1Tree.getCell(1, 2).outEdges.get(2);


    removedEdge1.removeThis();
    graph1Tree.edges.remove(removedEdge1);

    removedEdge2.removeThis();
    graph1Tree.edges.remove(removedEdge2);

    removedEdge3.removeThis();
    graph1Tree.edges.remove(removedEdge3);

    removedEdge4.removeThis();
    graph1Tree.edges.remove(removedEdge4);

    graph1.graphToTree();

    t.checkExpect(graph1, graph1Tree);




    Graph graph2Tree = new Graph(2, 2, new Random(2));

    t.checkExpect(graph2, graph2Tree);

    Edge removedEdge = graph2Tree.getCell(0, 0).outEdges.get(0);

    removedEdge.removeThis();
    graph2Tree.edges.remove(removedEdge);

    graph2.graphToTree();

    t.checkExpect(graph2, graph2Tree);

  }



  // Tests getCell
  void testGetCell(Tester t) {
    this.initData();
    t.checkExpect(this.kruskal4.getCell(3, 2), this.lCell);
    t.checkExpect(this.kruskal4.getCell(7, 8), this.rCell);
    t.checkExpect(this.kruskal4.getCell(1, 1), this.bCell);
    t.checkExpect(this.kruskal4.getCell(6, 7), this.openCell);


  }

  // Tests the MazeWorld constructors
  void testMazeWorldConstructors(Tester t) {
    t.checkExpect(world3, new MazeWorld(5, 5, new Random(3)));
    t.checkExpect(world4, new MazeWorld(10, 10, new Random(4)));
    t.checkExpect(world3.width, 5);
    t.checkExpect(world3.height, 5);
    t.checkExpect(world4.width, 10);
    t.checkExpect(world4.height, 10);
    t.checkExpect(world3.graph, kruskal3);
    t.checkExpect(world4.graph, kruskal4);

    Queue<Cell> breadthList3 = new ArrayDeque<>();
    Queue<Cell> breadthList4 = new ArrayDeque<>();
    Stack<Cell> depthList3 = new Stack<>();
    Stack<Cell> depthList4 = new Stack<>();

    breadthList3.add(world3.graph.getCell(0, 0));
    breadthList4.add(world4.graph.getCell(0, 0));
    depthList3.add(world3.graph.getCell(0, 0));
    depthList4.add(world4.graph.getCell(0, 0));


    t.checkExpect(world3.breadthWorkList, breadthList3);
    t.checkExpect(world4.breadthWorkList, breadthList4);

    t.checkExpect(world3.depthWorkList, depthList3);
    t.checkExpect(world4.depthWorkList, depthList4);

    t.checkExpect(world3.cameFromEdge, new HashMap<>());
    t.checkExpect(world4.cameFromEdge, new HashMap<>());

    t.checkExpect(world3.currCell, world3.graph.getCell(0, 0));
    t.checkExpect(world4.currCell, world4.graph.getCell(0, 0));

    t.checkExpect(world3.state, State.chooseState);
    t.checkExpect(world4.state, State.chooseState);

    t.checkExpect(world3.movesCount, 0);
    t.checkExpect(world4.movesCount, 0);

    t.checkExpect(world3.visitedPathVisible, true);
    t.checkExpect(world4.visitedPathVisible, true);

    MazeWorld randMaze = new MazeWorld(15, 20);

    t.checkExpect(randMaze.width, 15);
    t.checkExpect(randMaze.height, 20);
    t.checkExpect(randMaze.currCell, randMaze.graph.getCell(0, 0));
    t.checkExpect(randMaze.state, State.chooseState);
    t.checkExpect(randMaze.cameFromEdge, new HashMap<>());
    t.checkExpect(randMaze.movesCount, 0);
    t.checkExpect(randMaze.visitedPathVisible, true);

    MazeWorld randVerticalMaze = new MazeWorld(15, 20, 100000);

    t.checkExpect(randVerticalMaze.width, 15);
    t.checkExpect(randVerticalMaze.height, 20);
    t.checkExpect(randVerticalMaze.currCell, randVerticalMaze.graph.getCell(0, 0));
    t.checkExpect(randVerticalMaze.state, State.chooseState);
    t.checkExpect(randVerticalMaze.cameFromEdge, new HashMap<>());
    t.checkExpect(randVerticalMaze.movesCount, 0);
    t.checkExpect(randVerticalMaze.visitedPathVisible, true);

    MazeWorld randHorizontalMaze = new MazeWorld(15, 20, new Random(1), 0.000001);

    t.checkExpect(randHorizontalMaze.width, 15);
    t.checkExpect(randHorizontalMaze.height, 20);
    t.checkExpect(randHorizontalMaze.currCell, randHorizontalMaze.graph.getCell(0, 0));
    t.checkExpect(randHorizontalMaze.state, State.chooseState);
    t.checkExpect(randHorizontalMaze.cameFromEdge, new HashMap<>());
    t.checkExpect(randHorizontalMaze.movesCount, 0);
    t.checkExpect(randHorizontalMaze.visitedPathVisible, true);
    t.checkExpect(randHorizontalMaze.graph.getCell(0, 0).isConnected(
        randHorizontalMaze.graph.getCell(0, 1)), true);
  }

  // Tests makeScene 
  void testMakeScene(Tester t) {
    this.initData();

    WorldScene processScene = new WorldScene(50, 150);

    for (int row = 0; row < 5; row++) {
      for (int col = 0; col < 5; col++) {
        processScene.placeImageXY(this.graph3SearchInProg.getCell(
            row, col).drawCell(true), 5 + 10 * col, 5 + 10 * row);
      }
    }

    processScene.placeImageXY(new TextImage("moves: 0", 15, Color.BLACK), 25, 75);

    t.checkExpect(this.depthInProg.makeScene(), processScene);


    WorldScene chooseScene = new WorldScene(100, 200);

    for (int row = 0; row < 10; row++) {
      for (int col = 0; col < 10; col++) {
        chooseScene.placeImageXY(
            kruskal4.getCell(row, col).drawCell(true), 5 + 10 * col, 5 + 10 * row);
      }
    }

    chooseScene.placeImageXY(new TextImage("moves: 0", 15, Color.BLACK), 50 , 125);

    t.checkExpect(this.world4.makeScene(), chooseScene);


    WorldScene winScene = new WorldScene(50, 150);

    for (int row = 0; row < 5; row++) {
      for (int col = 0; col < 5; col++) {
        winScene.placeImageXY(
            this.graph3SolutionFin.getCell(row, col).drawCell(true), 5 + 10 * col, 5 + 10 * row);
      }
    }

    winScene.placeImageXY(new TextImage("moves: 0", 15, Color.BLACK), 25, 75);
    winScene.placeImageXY(new TextImage("91 % right", 15, Color.BLACK), 25, 125);

    t.checkExpect(this.finWin.makeScene(), winScene);

  }

  // Tests getRep method 
  void testGetRep(Tester t) {
    this.initData();
    t.checkExpect(new Utils().getRep(reps, midCell), this.midCell);
    t.checkExpect(new Utils().getRep(reps, leftCell), this.topCell);
    t.checkExpect(new Utils().getRep(reps, rightCell), this.topCell);
    t.checkExpect(new Utils().getRep(reps, bottomCell), this.midCell);
    t.checkExpect(new Utils().getRep(reps, topCell), this.topCell);
  }

  // tests the onKeyMethod
  void testOnKeyEvent(Tester t) {
    this.initData();
    MazeWorld chooseWorld1 = this.world3;
    this.initData();
    MazeWorld chooseWorld2 = this.world3;
    this.initData();
    MazeWorld chooseWorld3 = this.world3;
    
    chooseWorld1.onKeyEvent("b");
    chooseWorld2.onKeyEvent("d");
    chooseWorld3.onKeyEvent("p");
    chooseWorld1.onKeyEvent("p");
    chooseWorld2.onKeyEvent("b");
    chooseWorld3.onKeyEvent("d");
    
    t.checkExpect(chooseWorld1.state, State.breadthSearch);
    t.checkExpect(chooseWorld2.state, State.depthSearch);
    t.checkExpect(chooseWorld3.state, State.playerMove);
    t.checkExpect(chooseWorld3.currCell.state, CellState.currCell);
    
    this.depthInProg.onKeyEvent("c");
    
    Stack<Cell> depthWorkList2 = new Stack<>();
    Queue<Cell> breadthWorkList2 = new ArrayDeque<>();
    depthWorkList2.add(this.depthInProg.graph.getCell(0, 0));
    breadthWorkList2.add(this.depthInProg.graph.getCell(0, 0));
    
    t.checkExpect(this.depthInProg.state, State.chooseState);
    t.checkExpect(this.depthInProg.graph.getCell(0, 0).state, CellState.notVisited);
    t.checkExpect(this.depthInProg.graph.getCell(1, 0).state, CellState.notVisited);
    t.checkExpect(this.depthInProg.currCell, this.depthInProg.graph.getCell(0, 0));
    t.checkExpect(this.depthInProg.depthWorkList, depthWorkList2);
    t.checkExpect(this.depthInProg.breadthWorkList, breadthWorkList2);
    t.checkExpect(this.depthInProg.width, 5);
    t.checkExpect(this.depthInProg.height, 5);
    t.checkExpect(this.depthInProg.movesCount, 0);
    
    this.breadthInProg.onKeyEvent("r");
    
    Stack<Cell> depthWorkList3 = new Stack<>();
    Queue<Cell> breadthWorkList3 = new ArrayDeque<>();
    depthWorkList3.add(this.breadthInProg.graph.getCell(0, 0));
    breadthWorkList3.add(this.breadthInProg.graph.getCell(0, 0));
    
    t.checkExpect(this.breadthInProg.state, State.chooseState);
    t.checkExpect(this.breadthInProg.currCell, this.breadthInProg.graph.getCell(0, 0));
    t.checkExpect(this.breadthInProg.cameFromEdge, new HashMap<>());
    t.checkExpect(this.breadthInProg.depthWorkList, depthWorkList3);
    t.checkExpect(this.breadthInProg.breadthWorkList, breadthWorkList3);
    t.checkExpect(this.breadthInProg.width, 5);
    t.checkExpect(this.breadthInProg.height, 5);
    t.checkExpect(this.breadthInProg.movesCount, 0);
    
    t.checkExpect(this.animateWin.visitedPathVisible, true);
    this.animateWin.onKeyEvent("t");
    t.checkExpect(this.animateWin.visitedPathVisible, false);
    this.animateWin.onKeyEvent("t");
    t.checkExpect(this.animateWin.visitedPathVisible, true);
    
  }
 
  // tests the onTick() method
  void testOnTick(Tester t) {
    this.initData();

    MazeWorld breadthWorld = this.breadthInProg;
    MazeWorld depthWorld = this.depthInProg;
    MazeWorld winWorld = this.animateWin;

    this.initData();

    breadthWorld.onTick();
    depthWorld.onTick();
    winWorld.onTick();

    this.breadthInProg.searchBreadth();
    this.depthInProg.searchDepth();
    this.animateWin.animateWinPath();

    t.checkExpect(breadthWorld.movesCount, this.breadthInProg.movesCount);
    t.checkExpect(breadthWorld.breadthWorkList, this.breadthInProg.breadthWorkList);
    t.checkExpect(breadthWorld.graph, this.breadthInProg.graph);
    t.checkExpect(breadthWorld.state, this.breadthInProg.state);

    t.checkExpect(depthWorld.movesCount, this.depthInProg.movesCount);
    t.checkExpect(depthWorld.depthWorkList, this.depthInProg.depthWorkList);
    t.checkExpect(depthWorld.graph, this.depthInProg.graph);
    t.checkExpect(depthWorld.state, this.depthInProg.state);

    t.checkExpect(winWorld.movesCount, this.animateWin.movesCount);
    t.checkExpect(winWorld.depthWorkList, this.animateWin.depthWorkList);
    t.checkExpect(winWorld.graph, this.animateWin.graph);
    t.checkExpect(winWorld.state, this.animateWin.state);

  }


  // tests the searchDepth method  
  void testSearchDepth(Tester t) {
    this.initData();

    MazeWorld depth1 = this.depthInProg;
    MazeWorld depth2 = this.depthFin;

    this.initData();

    depth1.searchDepth();
    depth2.searchDepth();

    this.depthInProg.searchHelp(this.depthInProg.depthWorkList,
        this.depthInProg.depthWorkList.pop());
    this.depthFin.searchHelp(this.depthFin.depthWorkList,
        this.depthFin.depthWorkList.pop());

    t.checkExpect(depth1.movesCount, this.depthInProg.movesCount);
    t.checkExpect(depth1.depthWorkList, this.depthInProg.depthWorkList);
    t.checkExpect(depth1.graph, this.depthInProg.graph);
    t.checkExpect(depth1.state, this.depthInProg.state);

    t.checkExpect(depth2.movesCount, this.depthFin.movesCount);
    t.checkExpect(depth2.depthWorkList, this.depthFin.depthWorkList);
    t.checkExpect(depth2.graph, this.depthFin.graph);
    t.checkExpect(depth2.state, this.depthFin.state);
  }

  // tests the searchBreadth method 
  void testSearchBreadth(Tester t) {
    this.initData();

    MazeWorld breadth1 = this.breadthInProg;
    MazeWorld breadth2 = this.breadthFin;

    this.initData();

    breadth1.searchBreadth();
    breadth2.searchBreadth();

    this.breadthInProg.searchHelp(
        this.breadthInProg.breadthWorkList, this.breadthInProg.breadthWorkList.poll());
    this.breadthFin.searchHelp(
        this.breadthFin.breadthWorkList, this.breadthFin.breadthWorkList.poll());

    t.checkExpect(breadth1.movesCount, this.breadthInProg.movesCount);
    t.checkExpect(breadth1.breadthWorkList, this.breadthInProg.breadthWorkList);
    t.checkExpect(breadth1.graph, this.breadthInProg.graph);
    t.checkExpect(breadth1.state, this.breadthInProg.state);

    t.checkExpect(breadth2.movesCount, this.breadthFin.movesCount);
    t.checkExpect(breadth2.breadthWorkList, this.breadthFin.breadthWorkList);
    t.checkExpect(breadth2.graph, this.breadthFin.graph);
    t.checkExpect(breadth2.state, this.breadthFin.state);

  }

  // tests the searchHelp(ArrayList  workList, Cell c) method 
  void testSearchHelp(Tester t) {
    this.initData();

    this.breadthInProg.searchHelp(
        this.breadthInProg.breadthWorkList, this.breadthInProg.breadthWorkList.poll());

    t.checkExpect(this.breadthInProg.movesCount, 1);
    t.checkExpect(this.breadthInProg.breadthWorkList.contains(
        this.breadthInProg.graph.getCell(1, 1)), false);
    t.checkExpect(this.breadthInProg.breadthWorkList.contains(
        this.breadthInProg.graph.getCell(2, 1)), true);
    t.checkExpect(this.breadthInProg.graph.getCell(1, 1).state,
        CellState.visited);
    t.checkExpect(
        this.breadthInProg.cameFromEdge.keySet().contains(
            this.breadthInProg.graph.getCell(2, 1)), true);

    this.breadthFin.searchHelp(this.breadthFin.breadthWorkList,
        this.breadthFin.breadthWorkList.poll());

    t.checkExpect(this.breadthFin.movesCount, 1);
    t.checkExpect(this.breadthFin.breadthWorkList.contains(
        this.breadthFin.graph.getCell(4, 4)), false);
    t.checkExpect(
        this.breadthFin.graph.getCell(4, 4).state, CellState.visited);
    t.checkExpect(
        this.breadthFin.cameFromEdge.keySet().contains(
            this.breadthFin.graph.getCell(4, 4)), true);
    t.checkExpect(this.breadthFin.state, State.winState);
    t.checkExpect(this.breadthFin.currCell, this.breadthFin.graph.getCell(4, 4));




    this.depthInProg.searchHelp(this.depthInProg.depthWorkList, 
        this.depthInProg.depthWorkList.pop());

    t.checkExpect(this.depthInProg.movesCount, 1);
    t.checkExpect(this.depthInProg.depthWorkList.contains(
        this.depthInProg.graph.getCell(3, 0)), false);
    t.checkExpect(this.depthInProg.depthWorkList.contains(
        this.depthInProg.graph.getCell(4, 0)), true);
    t.checkExpect(this.depthInProg.graph.getCell(3, 0).state, CellState.visited);
    t.checkExpect(this.depthInProg.cameFromEdge.keySet().contains(
        this.depthInProg.graph.getCell(4, 0)), true);

    this.depthFin.searchHelp(this.depthFin.depthWorkList, 
        this.depthFin.depthWorkList.pop());

    t.checkExpect(this.depthFin.movesCount, 1);
    t.checkExpect(this.depthFin.depthWorkList.contains(
        this.depthFin.graph.getCell(4, 4)), false);
    t.checkExpect(this.depthFin.graph.getCell(4, 4).state,
        CellState.visited);
    t.checkExpect(this.breadthFin.cameFromEdge.keySet().contains(
        this.breadthFin.graph.getCell(4, 4)), true);
    t.checkExpect(this.depthFin.state, State.winState);
    t.checkExpect(this.depthFin.currCell, this.depthFin.graph.getCell(4, 4));

  }


  // tests the reset method 
  void testReset(Tester t) {
    this.initData();

    this.graph3SearchInProg.reset();
    this.graph3SolutionProg.reset();

    t.checkExpect(this.graph3SearchInProg.getCell(0, 0).state, CellState.notVisited);
    t.checkExpect(this.graph3SearchInProg.getCell(2, 0).state, CellState.notVisited);
    t.checkExpect(this.graph3SolutionProg.getCell(0, 1).state, CellState.notVisited);
    t.checkExpect(this.graph3SolutionProg.getCell(4, 4).state, CellState.notVisited);

  }

  // tests the countOnPath() Method 
  void testCountOnPath(Tester t) {
    this.initData();

    t.checkExpect(this.graph3SearchInProg.countOnPath(), 0);
    t.checkExpect(this.graph3SearchEnd.countOnPath(), 0);
    t.checkExpect(this.graph3SolutionProg.countOnPath(), 5);
    t.checkExpect(this.graph3SolutionFin.countOnPath(), 11);
  }

  // tests the countHasVisited() method
  void testCountHasVisited(Tester t) {
    this.initData();

    t.checkExpect(this.graph3SearchInProg.countHasVisited(), 3);
    t.checkExpect(this.graph3SearchEnd.countHasVisited(), 11);
    t.checkExpect(this.graph3SolutionProg.countHasVisited(), 12);
    t.checkExpect(this.graph3SolutionFin.countHasVisited(), 12);
  }

  
  // tests the cellColor(boolean b) method 
  void testCellColor(Tester t) {
    this.initData();
    t.checkExpect(cellNotVisited.cellColor(false), Color.GRAY);
    t.checkExpect(cellNotVisited.cellColor(true), Color.gray);
    t.checkExpect(cellVisited.cellColor(true), Color.cyan);
    t.checkExpect(cellVisited.cellColor(false), Color.gray);
    t.checkExpect(cellOnPath.cellColor(true), Color.BLUE);
    t.checkExpect(cellOnPath.cellColor(false), Color.blue);
    t.checkExpect(cellCurrCell.cellColor(false), Color.red);
    t.checkExpect(cellCurrCell.cellColor(true), Color.red);
  }

  // tests the hasVisited() method
  void testHasVisited(Tester t) {
    this.initData();
    t.checkExpect(cellNotVisited.hasVisited(), false);
    t.checkExpect(cellVisited.hasVisited(), true);
    t.checkExpect(cellOnPath.hasVisited(), false);
    t.checkExpect(cellCurrCell.hasVisited(), false);
  }

  // tests the isOnPath() method
  void testIsOnPath(Tester t) {
    this.initData();
    t.checkExpect(cellNotVisited.isOnPath(), false);
    t.checkExpect(cellVisited.isOnPath(), false);
    t.checkExpect(cellOnPath.isOnPath(), true);
    t.checkExpect(cellCurrCell.isOnPath(), false);
  }


  // tests the cellAtPos(int r, int col) method 
  void testCellAtPos(Tester t) {
    this.initData();
    t.checkExpect(this.cellNotVisited.cellAtPos(3, 2), true);
    t.checkExpect(this.lCell.cellAtPos(5, 2), false);
    t.checkExpect(this.cellVisited.cellAtPos(7, 8), true);
    t.checkExpect(this.bCell.cellAtPos(10, 9), false);

  }

  // tests the setVisited() method 
  void testSetVisited(Tester t) {
    this.initData();
    this.cellVisited.setVisited();
    t.checkExpect(cellVisited.state, CellState.visited);
    this.cellNotVisited.setVisited();
    t.checkExpect(cellNotVisited.state, CellState.visited);
    this.cellOnPath.setVisited();
    t.checkExpect(cellOnPath.state, CellState.visited);
    this.cellCurrCell.setVisited();
    t.checkExpect(cellCurrCell.state, CellState.visited);

  }

  // tests the setNotVisited() method 
  void testSetNotVisited(Tester t) {
    this.initData();
    this.cellVisited.setNotVisited();
    t.checkExpect(cellVisited.state, CellState.notVisited);

    this.cellOnPath.setNotVisited();
    t.checkExpect(cellOnPath.state, CellState.notVisited);

    this.cellCurrCell.setNotVisited();
    t.checkExpect(cellCurrCell.state, CellState.notVisited);
  }

  // tests the setCurrCell() method 
  void testSetCurrCell(Tester t) {
    this.initData();
    this.cellVisited.setCurrCell();
    t.checkExpect(cellVisited.state, CellState.currCell);
    this.cellOnPath.setCurrCell();
    t.checkExpect(cellVisited.state, CellState.currCell);
    this.cellNotVisited.setCurrCell();
    t.checkExpect(cellNotVisited.state, CellState.currCell);
    this.cellCurrCell.setCurrCell();
    t.checkExpect(cellCurrCell.state, CellState.currCell);


  }

  // tests the setOnPath() method
  void testSetOnPath(Tester t) {
    this.initData();
    this.cellVisited.setOnPath();
    t.checkExpect(cellVisited.state, CellState.onPath);
    this.cellNotVisited.setOnPath();
    t.checkExpect(cellNotVisited.state, CellState.onPath);
    this.cellCurrCell.setOnPath();
    t.checkExpect(cellCurrCell.state, CellState.onPath);
    this.cellOnPath.setOnPath();
    t.checkExpect(cellOnPath.state, CellState.onPath);

  }

  // tests the moveCell(Graph graph, int r, int c) method 
  void testMoveCell(Tester t) {
    this.initData();
    t.checkExpect(this.lCell.moveCell(kruskal4, 2, 3), kruskal4.getCell(5, 5));
    //checking out of bounds 
    t.checkExpect(this.lCell.moveCell(kruskal4, 10, 30), null);
    t.checkExpect(this.rCell.moveCell(kruskal4, 1, 2), kruskal4.getCell(8, 10));

  }

  // tests the isConnected(Cell c) method
  void testIsConnected(Tester t) {
    this.initData();
    t.checkExpect(this.midCell.isConnected(topCell), true);
    t.checkExpect(this.midCell.isConnected(topCell), true);
    t.checkExpect(this.topCell.isConnected(bottomCell), false);
  }
  

  // tests the connectedEdge(Cell c) method
  void testConnectedEdge(Tester t) {
    this.initData();

    t.checkExpect(this.midCell.connectedEdge(this.leftCell), this.edgeLeft);
    t.checkExpect(this.midCell.connectedEdge(bottomCell), this.edgeBottom);
    t.checkExpect(this.rightCell.connectedEdge(bottomCell), null);
  }

  // tests the getOtherCell(Cell c) method
  void testGetOtherCell(Tester t) {
    this.initData();
    t.checkExpect(this.edgeLeft.getOtherCell(midCell), this.leftCell);
    t.checkExpect(this.edgeBottom.getOtherCell(this.bottomCell), this.midCell);
    t.checkException(new RuntimeException("cell not next to edge"), this.edgeLeft, 
        "getOtherCell", bottomCell);

  }

  // tests the forEachEdge(Consumer<Edge> func) method 
  void testForEachEdge(Tester t) {
    this.initData();
    this.leftCell.forEachEdge((edge -> {
      edge.removeThis(); } ));
    t.checkExpect(this.leftCell.outEdges.contains(edgeLeft), false);
    t.checkExpect(this.midCell.outEdges.contains(edgeLeft), false);
    this.initData();
    t.checkExpect(this.leftCell.outEdges.contains(edgeLeft), true);
    t.checkExpect(this.leftCell.outEdges.contains(edgeLeft), true);

  }

  // tests the moveInDirection(Cell c) method 
  void testMoveInDirection(Tester t) {
    this.initData();
    this.midCell.setCurrCell();
    this.world4.moveInDirection(leftCell);
    t.checkExpect(this.leftCell.state, CellState.currCell);


    this.initData();
    this.world4.moveInDirection(midCell);
    t.checkExpect(this.midCell.state, CellState.currCell);

    // test null-- not connected
    this.initData();
    this.world4.moveInDirection(null);
    t.checkExpect(this.bottomCell.state, CellState.currCell);

  }

  // tests the move() method 
  void testMove(Tester t) {
    this.initData();
    this.midCell.setCurrCell();
    this.world4.move("up");
    t.checkExpect(this.topCell.state, CellState.currCell);

    this.initData();
    this.midCell.setCurrCell();
    this.world4.move("left");
    t.checkExpect(this.leftCell.state, CellState.currCell);

    this.initData();
    this.midCell.setCurrCell();
    this.world4.move("right");
    t.checkExpect(this.rightCell.state, CellState.currCell);

  }
  
  // tests the animateWithPath() method
  void testAnimateWinPath(Tester t) {

    this.initData();
    world3.animateWinPath();
    t.checkExpect(world3.currCell, world3.graph.getCell(0,0));
    world4.animateWinPath();
    t.checkExpect(world4.currCell, world4.graph.getCell(0,0));

  }

}


