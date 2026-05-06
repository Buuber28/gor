public  class Gor4ModelData {

    public final String aaOrder;
    public final long[][][][][][] counts6D;
    public final long[][][][] counts4D;


    public Gor4ModelData(String aaOrder) {
        this.aaOrder = aaOrder;
        this.counts6D = new long[3][20][20][17][20][17];
        this.counts4D = new long[3][20][20][17];

    }
}