class cgpWrapper {
public:
    void harmonic_runCGP();
    static void my_runCGP();
private:
    static double fitness(struct parameters*, struct chromosome*, struct dataSet*);
};