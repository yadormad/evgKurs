package umlKurs;

public class Lambda {



    private double f(double x, double l, double H) {
        return f1(x) - f2(x, l, H);
    }

    private double f1(double x) {
        return Math.tan(x);
    }

    private double f2(double x, double l, double H) {
        return H*l/x;
    }

    private boolean isEquals(double y1, double y2, double e){
        if(Math.abs(y1 - y2) < e)
            return true;
        else
            return false;
    }

    /*private double halfDivideMethod(double a, double b, double ) {

    }
*/
}
