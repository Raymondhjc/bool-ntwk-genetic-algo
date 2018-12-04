import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

public class Program {

    final static int GAMA = 256;
    final static int H = 3;
    final static int ITERATIONS = 1000;
    final static double MUTATION_PROB = 0.01;
    final static double MIN_MUTUAL_INFO = 0.05;

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

        double[] dynamicAccuracy = new double[names.size()];
        for(int nameIndex = 0; nameIndex < names.size(); ++nameIndex) {
            String st;
            try {
                br = new BufferedReader(new FileReader("files/" + names.get(nameIndex)));

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

                    // TODO one way to k mean
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

                    // calculate mutual information
                    double[][] MI = mi(b);

                    int[][] network = new int[N][];
                    // start genetic algorithm
                    double curDA = 0;
                    for(int i = 0; i < N; ++i) {
                        network[i] = genetic(b, POPULATION, MAX_REG_NUM, MI[i], i);
//                        System.out.println("for gene G" + i);
//                        for(int g : network[i]) {
//                            System.out.print(g + " ");
//                        }
//                        System.out.print("\n");
                        curDA += consistency(network[i], b, i);
                    }
                    dynamicAccuracy[nameIndex] = Math.max(curDA, dynamicAccuracy[nameIndex]);
                }

            } catch(Exception e) {
                e.printStackTrace();
            }
        }
        for(int i = 0; i < names.size(); ++i) {
            System.out.println(names.get(i).split("-")[1]);
            System.out.println("Dynamics consistency is : ---- " + dynamicAccuracy[i]);
        }
    }

    private static int[] genetic(boolean[][] data, int population, int maxReg, double[] mi, int v) {

        // if the gene v is a self-regulatory gene
        int count = 0, self = 0;
        for(int i = 0; i < mi.length; ++i) {
            if(mi[i] >= MIN_MUTUAL_INFO) {
                ++count;
                self = i;
            }
        }
        if(count <= 1 && self == v) {
            return new int[]{v};
        }

        // initial population
        int[][] pop = new int[population][];
        for(int i = 0; i < population; ++i) {
            int[] chromo = genInitialChromo(maxReg, data[0].length);
            // TODO don't do MI now
            if(exist(pop, chromo) || !isChomoValid(chromo, maxReg, mi)) {
                --i;
            } else {
                pop[i] = chromo;
            }
        }

        // initial population fitness value
        double[] fitness = new double[population];
        for(int i = 0; i < population; ++i) {
            fitness[i] = fitness(pop[i], data, v);
        }

        // GA iteration
        for(int it = 0; it < ITERATIONS; ++it) {

            // selection
            // calculate adjusted fitness value
            double FMAX = Double.MIN_VALUE, FMIN = Double.MAX_VALUE;
            for(double f : fitness) {
                FMAX = Math.max(FMAX, f);
                FMIN = Math.min(FMIN, f);
            }

            // TODO convergence
            if(FMAX == FMIN) {
                // return the one with best fitness
                int index = 0;
                double max = Double.MIN_VALUE;
                for(int i = 0; i < population; ++i) {
                    if(fitness[i] > max) {
                        max = fitness[i];
                        index = i;
                    }
                }
//                System.out.println("iterations:" + it);
                return pop[index];
            }

            double[] adjFitness = new double[population];

            // do roulette wheel selection
            double sum = 0;
            for(int i = 0; i < population; ++i) {
                double a = H / (FMAX - FMIN);
                adjFitness[i] = a * fitness[i] + 1 - a * FMIN;
                sum += adjFitness[i];
            }

            // parents p1 and p2, index of chromosomes (pop[])
            int p1 = 0;
            Random rand = new Random();
            double r = rand.nextDouble();
            for(int i = 0, acc = 0; i < population; ++i) {
                if(r <= (acc += adjFitness[i]) / sum) {
                    p1 = i;
                    break;
                }
            }
            int p2 = p1;
            while(p2 == p1) {
                r = rand.nextDouble();
                for(int i = 0, acc = 0; i < population; ++i) {
                    if(r <= (acc += adjFitness[i]) / sum) {
                        p2 = i;
                        break;
                    }
                }
            }

            // crossover and mutation
            int n = data[0].length;
            boolean[] parent1 = new boolean[n];
            boolean[] parent2 = new boolean[n];
            for(int i = 0, i1 = 0, i2 = 0; i < n; ++i) {
                if(i1 < pop[p1].length && i == pop[p1][i1]) {
                    parent1[i] = true;
                    ++i1;
                }
                if(i2 < pop[p2].length && i == pop[p2][i2]) {
                    parent2[i] = true;
                    ++i2;
                }
            }

            boolean[] offSpring1;
            boolean[] offSpring2;
            int n1 = 0, n2 = 0;
            offSpring1 = new boolean[n];
            offSpring2 = new boolean[n];
            for(int i = 0; i < n; ++i) {
                if(parent1[i] == parent2[i]) {
                    offSpring1[i] = parent1[i];
                    offSpring2[i] = parent1[i];
                } else {
                    offSpring1[i] = rand.nextInt(2) == 1;
                    offSpring2[i] = rand.nextInt(2) == 1;
                }
                // mutation
                offSpring1[i] = rand.nextDouble() <= MUTATION_PROB != offSpring1[i];
                offSpring2[i] = rand.nextDouble() <= MUTATION_PROB != offSpring2[i];

                if(offSpring1[i]) {
                    ++n1;
                }
                if(offSpring2[i]) {
                    ++n2;
                }
            }

            int[] c1 = new int[n1];
            int[] c2 = new int[n2];
            for(int i = 0, i1 = 0, i2 = 0; i < n; ++i) {
                if(offSpring1[i]) {
                    c1[i1++] = i;
                }
                if(offSpring2[i]) {
                    c2[i2++] = i;
                }
            }

            // replacement
            double offSpring1Fitness = fitness(c1, data, v);
            double offSpring2Fitness = fitness(c2, data, v);
            if(isChromoValid(offSpring1, maxReg, mi) && offSpring1Fitness > fitness[p1]) {
                pop[p1] = c1;
                fitness[p1] = offSpring1Fitness;
            }
            if(isChromoValid(offSpring2, maxReg, mi) && offSpring2Fitness > fitness[p2]) {
                pop[p2] = c2;
                fitness[p2] = offSpring2Fitness;
            }
        }

        // return the one with best fitness
        int index = 0;
        double max = Double.MIN_VALUE;
        for(int i = 0; i < population; ++i) {
            if(fitness[i] > max) {
                max = fitness[i];
                index = i;
            }
        }

        return pop[index];
    }

    private static int[] genInitialChromo(int maxReg, int len) {
        // at least 1 and at most maxReg regulatory genes
        int num = 0;
        while(num == 0) {
            num = new Random().nextInt(maxReg) + 1;
        }
        int[] chromosome = new int[num];
        Set<Integer> set = new HashSet<>();
        for(int i = 0; i < num; ++i) {
            Random r = new Random();
            int rand = r.nextInt(len);
            // can be self regulatory
            while(set.contains(rand)) {
                rand = r.nextInt(len);
            }
            set.add(rand);
            chromosome[i] = rand;
        }
        Arrays.sort(chromosome);
        return chromosome;
    }

    private static boolean exist(int[][] pop, int[] tmp) {
        for(int i = 0; i < pop.length && pop[i] != null; ++i) {
            boolean exist = true;
            if(pop[i].length == tmp.length) {
                for(int j = 0; j < pop[i].length; ++j) {
                    if(tmp[j] != pop[i][j]) {
                        exist = false;
                        break;
                    }
                }
            } else {
                exist = false;
            }
            if(exist) return true;
        }
        return false;
    }

    // is the candidate chromosome valid, check max regulatory gene number and mutual information
    private static boolean isChomoValid(int[] chromo, int maxReg, double[] mi) {
//        if(chromo.length > maxReg) {
//            return false;
//        }
//        for(int v : chromo) {
//            if(mi[v] < MIN_MUTUAL_INFO) {
//                return false;
//            }
//        }
        return true;
    }
    private static boolean isChromoValid(boolean[] chromo, int maxReg, double[] mi) {
//        int cnt = 0;
//        for(boolean b : chromo) {
//            if(b) {
//                ++cnt;
//            }
//        }
//        if(cnt > maxReg) {
//            return false;
//        }
//        for(int i = 0; i < chromo.length; ++i) {
//            if(chromo[i] && mi[i] < MIN_MUTUAL_INFO) {
//                return false;
//            }
//        }
        return true;
    }

    private static double fitness(int[] chromo, boolean[][] data, int v) {
        double C = consistency(chromo, data, v);
        return 1 / ((1 - C) * GAMA + chromo.length);
    }

    private static double consistency(int[] chromo, boolean[][] data, int v) {
        Map<Integer, Integer> map= new HashMap<>();
        for(int i = 0; i < data.length - 1; ++i) {
            int b = 0;
            for(int c : chromo) {
                b <<= 1;
                if(data[i][c]) {
                    b += 1;
                }
            }
            int prediction = data[i + 1][v] ? 1 : -1;
            map.put(b, map.getOrDefault(b, 0) + prediction);
        }

        //calculate consistency, from time t = 2
        int count = 0;
        for(int i = 0; i < data.length - 1; ++i) {
            int b = 0;
            for(int c : chromo) {
                b <<= 1;
                if(data[i][c]) {
                    b += 1;
                }
            }
            if(map.get(b) > 0 && data[i + 1][v] || map.get(b) < 0 && !data[i + 1][v]) {
                ++count;
            }
        }
        return (double)count / (data.length - 1);
    }

    // calculate mutual information
    private static double[][] mi(boolean[][] data) {
        int m = data.length, n = data[0].length;
        int[][] p = new int[n][2];
        // conditional probability : 0|0, 0|1, 1|0, 1|1
        int[][][] cond = new int[n][n][4];
        for(boolean[] b : data) {
            for(int j = 0; j < n; ++j) {
                if(b[j]) {
                    ++p[j][1];
                    for(int k = 0; k < n; ++k) {
                        if(b[k]) {
                            // 11
                            ++cond[j][k][3];
                        } else {
                            // 10
                            ++cond[j][k][2];
                        }
                    }
                } else {
                    ++p[j][0];
                    for(int k = 0; k < n; ++k) {
                        if(b[k]) {
                            // 01
                            ++cond[j][k][1];
                        } else {
                            // 00
                            ++cond[j][k][0];
                        }
                    }
                }
            }
        }
        double[][] mi = new double[n][n];
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                for(int k = 0; k < 4; ++k) {
                    if(cond[i][j][k] != 0) {
                        mi[i][j] += (double)cond[i][j][k] / m
                                * Math.log((double)cond[i][j][k] / p[i][0] / p[i][1] * m);
                    }
                }
            }
        }
        return mi;
    }

}
