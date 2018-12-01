import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

public class Program {

    final static int GAMA = 256;
    final static int H = 3;
    final static int ITERATIONS = 1000;
    // mutation probability, 1/PROB
    final static int MUTATION_PROB = 100;

    public static void main(String args[]) {
        File file = new File("config.cfg");
        BufferedReader br;
        // file names
        List<String> names = new ArrayList<>();
        try {
            String st;
            br = new BufferedReader(new FileReader(file));
            while((st = br.readLine()) != null) {
                names.add(st);
            }
        } catch(Exception e) {
            e.printStackTrace();
        }

        for(String name : names) {
            String st;
            try {
                br = new BufferedReader(new FileReader("files/" + name));

                // read head
                st = br.readLine();
                // N = genes count
                final int N = st.trim().split("\\s+").length - 1;
                // population size
                final int POPULATION = N + 10;
                // Maximum number of regulatory genes
                final int MAX_REG_NUM = (int)Math.floor(N * 0.6);

                while((st = br.readLine()) != null) {

                    // read expression volume
                    List<double[]> exprs = new ArrayList<>();
                    double means[] = new double[N];
                    while(st != null && !st.equals("")) {
                        String[] s = st.trim().split("\\s+");
                        double[] expr = new double[N];
                        for(int i = 0; i < N; ++i) {
                            expr[i] = Double.valueOf(s[i + 1]);
                            means[i] += expr[i];
                        }
                        exprs.add(expr);
                        st = br.readLine();
                    }

                    // get average
                    for(int i = 0; i < means.length; ++i) {
                        means[i] /= exprs.size();
                    }

                    // to boolean data set
                    boolean[][] b = new boolean[exprs.size()][N];
                    for(int i = 0; i < exprs.size(); ++i) {
                        double[] e = exprs.get(i);
                        for(int j = 0; j < e.length; ++j) {
                            b[i][j] = e[j] > means[j];
                        }
                    }

                    // start genetic algorithm
                    genetic(b);
                }

            } catch(Exception e) {
                e.printStackTrace();
            }
        }
    }

    private static void genetic(boolean[][] dataset) {

    }

}
