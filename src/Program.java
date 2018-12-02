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
                    for(int i = 0; i < N; ++i) {
                        genetic(b, POPULATION, MAX_REG_NUM, i);
                    }
                }

            } catch(Exception e) {
                e.printStackTrace();
            }
        }
    }

    private static void genetic(boolean[][] data, int population, int maxReg, int v) {
        // initial population
        int[][] pop = new int[population][];
        for(int i = 0; i < population; ++i) {
            int[] tmp = genChromo(maxReg, data.length, v);
            if(exist(pop, tmp)) {
                --i;
            } else {
                pop[i] = tmp;
            }
        }

        // initial population fitness value
        int[] fitness = new int[population];
        for(int i = 0; i < population; ++i) {
            fitness[i] = fitness(pop[i], data, v);
        }

        // GA
        for(int it = 0; it < ITERATIONS; ++it) {

            // selection

        }

    }

    private static int[] genChromo(int maxReg, int len, int v) {
        int num = new Random().nextInt(maxReg) + 1;
        int[] chromosome = new int[num];
        Set<Integer> set = new HashSet<>();
        for(int i = 0; i < num; ++i) {
            Random r = new Random();
            int rand = r.nextInt(len);
            while(set.contains(rand) || rand == v) {
                rand = r.nextInt(len);
            }
            set.add(rand);
            chromosome[i] = rand;
        }
        Arrays.sort(chromosome);
        return chromosome;
    }

    private static boolean exist(int[][] pop, int[] tmp) {
        for(int[] p : pop) {
            boolean f = false;
            for(int i = 0; i < pop.length; ++i) {
                if(tmp[i] != p[i]) {
                    f = true;
                    break;
                }
            }
            if(!f) {
                return false;
            }
        }
        return true;
    }

    private static int fitness(int[] chromo, boolean[][] data, int v) {
        for(int i = 1; i < data.length; ++i) {
            BitSet b = new BitSet(chromo.length);
        }
    }

}
