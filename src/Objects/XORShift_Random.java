package Objects;

public class XORShift_Random {
    private long seed;
    private long cnt;

    public XORShift_Random() {
        this.seed = System.nanoTime() | 1;
        cnt = seed;
    }


    public int getInt(int bound) {
        seed ^= (seed << 21);
        seed ^= (seed >>> 35);
        seed ^= (seed << 4);
        cnt += 123456789123456789L;
        int result = (int) ((seed + cnt) % bound);
        return (result < 0) ? -result : result;
    }
}
