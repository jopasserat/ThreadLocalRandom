import fr.isima.random.ThreadLocalMRG32k3a;

/**
 *
 * @author jopasserat
 */
public class Main implements Runnable {

    @Override
    public void run() {
        
        for (int i = 0; i < 5; ++i) {
            System.out.println("Thread[" +
                            ThreadLocalMRG32k3a.threadId.get() + 
                            "] {" + i + "} = " + ThreadLocalMRG32k3a.current().next());
        }
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Thread tArray[] = new Thread[3];
        
        for (Thread t : tArray) {
            t = new Thread(new Main());
            t.start();
        }
    }
}
