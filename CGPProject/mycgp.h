class cgpWrapper {
public:
    static void my_runCGP();
private:
    static double fitness(struct parameters*, struct chromosome*, struct dataSet*);
};