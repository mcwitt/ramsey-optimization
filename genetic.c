/* Genetic algorithm based on Goldberg's "Simple Genetic Algorithm" */
#define MAXPOP      100
#define MAXSTRING   30


typedef struct
{
    int chrom[MAXSTRING];
    double fitness;
} individual;

/* find the leftmost insertion point for x in a sorted array */
int bisect(double a[], double x, int l, int r)
{
    int mid;

    while (r-l > 1)
    {
        mid = (l+r)/2;
        if (a[mid] < x) l = mid;
        else r = mid;
    }

    return (x < a[l]) ? l : r;
}

void select(individual pop[], int popsize, double sumfitness)
{
}

/* do one-point crossover at position k */
void crossover(int i1[], int i2[], int o1[], int o2[], int k)
{
}

/* mutate string s with bit-flipping probability p_flip */
void mutate(int s[], double p_flip)
{
}
