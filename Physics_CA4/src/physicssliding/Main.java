
package physicssliding;
import java.util.Vector;

public class Main
{
    public static void main(String[] args)
    {
        //User must be able to set normal values
        
        Vector<Double> fNet=new Vector(3);
        fNet.add(0.0);
        fNet.add(0.0);
        fNet.add(0.0);
        
        int stepNum=3;
        double startTime=10.0;
        double finalTime=10.3;
        
        double stepSize=((finalTime-startTime)/stepNum);
        
        //____________________________________
        // variables for the plane
        //____________________________________
        Vector<Double> normal = new Vector(3);
        normal.add(2.0);
        normal.add(-4.0);
        normal.add(3.0);
        double mewStatic = 0.7;
        double mewKinetic = 0.4;
        //____________________________________
        //variables for the object
        //____________________________________
        Object object = new Object(0.8);
        Vector<Double>position=new Vector(3);
        Vector<Double>velocity=new Vector(3);
        Vector<Double>acceleration=new Vector(3);
        position.add(0.0);
        position.add(0.0);
        position.add(0.0);
        velocity.add(0.0);
        velocity.add(0.0);
        velocity.add(0.0);            
        //____________________________________
        
        //Force applied by the plane onto the object
        Vector nDirection = findNormalDirection(normal);
        
        //force of gravity
        Vector forceOfGravity = calculateGravityForce(object.getMass());
        
        //ammount of gravity in the direction of the normal force
        Vector Fgn = calculateFgn(forceOfGravity, nDirection);
        
        //ammount of gravity in the direction of the plane
        Vector Fgp = subtractVectors(forceOfGravity, Fgn);
        
        //how strong the force of gravity is in the direction of the plane
        double fgpLength = Math.sqrt(dotProduct(Fgp, Fgp));
        
        //the normal force acting on the object
        Vector Fn = calculateFn(Fgn);
        
        //Magnitude of the normal force
        double fnLength = Math.sqrt(dotProduct(Fn, Fn));
        
        
        System.out.println("________________________________________");
        System.out.println("N Direction: " + nDirection);
        System.out.println("Force of gravity: " + forceOfGravity);
        System.out.println("Fgn: " + Fgn);
        System.out.println("Fgp: " + Fgp);
        System.out.println("Fgp length: " + fgpLength + " N");
        System.out.println("Fn: " + Fn);
        System.out.println("fnLength: " + fnLength);
        System.out.println("________________________________________");
        System.out.println("");
        
        //See if the force of gravity in the direction of the plane is enough to move the object
        boolean moving=checkIfMoving(mewStatic, fnLength, fgpLength);
        if(moving)
        {
            Vector FkfDirection=findFrictionDirection(Fgp,fgpLength);
            double FkfLength=mewKinetic*fnLength;
            Vector kineticFriction=calculateFrictionForce(FkfLength, FkfDirection);
            fNet=addVectors(kineticFriction, Fgp);
            System.out.println("________________________________________");
            System.out.println("Ff Direction: "+FkfDirection);
            System.out.println("Ff Length: "+FkfLength);
            System.out.println("Ff: "+kineticFriction);
            System.out.println("Fnet: "+fNet);
            System.out.println("________________________________________");
            acceleration=calculateAcceleration(object, fNet);
            System.out.println("Acceleration: "+acceleration);            
            System.out.println("________________________________________");
            eulersMethod(object, position, velocity, fNet, stepNum, stepSize, startTime,  finalTime);
    
        }
        else
        {
            
        }            
    }

    public static double dotProduct(Vector<Double> V1, Vector<Double> V2)
    {
        double i;
        double j;
        double k;

        i = V1.get(0) * V2.get(0);
        j = V1.get(1) * V2.get(1);
        k = V1.get(2) * V2.get(2);

        //dotProduct = Math.sqrt(i + j + k);
        double dotProduct = i + j + k;

        return dotProduct;
    }

    public static Vector findNormalDirection(Vector<Double> normal)
    {
        Vector<Double> vector = new Vector<>();

        double dotProduct = dotProduct(normal, normal);

        vector.add(normal.get(0) * (1/Math.sqrt(dotProduct)));
        vector.add(normal.get(1) * (1/Math.sqrt(dotProduct)));
        vector.add(normal.get(2) * (1/Math.sqrt(dotProduct)));

        return vector;
    }

    public static Vector<Double> calculateGravityForce(double mass)
    {
        //Fg = -mgk
        double gravity = 9.81;
        double i;
        double j;
        double k;

        mass *= -1;

        Vector<Double> gravityForce = new Vector<>(3);

        i = 0;
        j = 0;
        k = mass * gravity;

        gravityForce.add(0, i);
        gravityForce.add(1, j);
        gravityForce.add(2, k);

        return gravityForce;
    }

    public static Vector<Double> calculateFgn(Vector<Double> Fg, Vector<Double> normalDirection)
    {
     //Fgn variable name goes against naming conventions - worth it for readability?
     Vector<Double> Fgn = new Vector(3);

     double dotProduct = dotProduct(Fg, normalDirection);

     Fgn.add(dotProduct * normalDirection.get(0));
     Fgn.add(dotProduct * normalDirection.get(1));
     Fgn.add(dotProduct * normalDirection.get(2));

     return Fgn;
 }

    public static Vector<Double> subtractVectors(Vector<Double> velocity, Vector<Double> wind)
    {
        Vector<Double> answer = new Vector(3);

        double velocityI = velocity.get(0);
        double velocityJ = velocity.get(1);
        double velocityK = velocity.get(2);

        double windI = wind.get(0);
        double windJ = wind.get(1);
        double windK = wind.get(2);

        answer.add(0, velocityI - windI);
        answer.add(1, velocityJ - windJ);
        answer.add(2, velocityK - windK);

        return answer;
    }

    public static Vector<Double> addVectors(Vector<Double>F1, Vector<Double>F2)
    {
      Vector<Double> answer = new Vector(3);

    double F1I = F1.get(0);
    double F1J = F1.get(1);
    double F1K = F1.get(2);

    double F2I = F2.get(0);
    double F2J = F2.get(1);
    double F2K = F2.get(2);

    answer.add(0, F1I + F2I);
    answer.add(1, F1J + F2J);
    answer.add(2, F1K + F2K);

        return answer;
  }

    public static Vector<Double> calculateFn(Vector<Double> Fgn)
    {
      Vector<Double> Fn = new Vector(3);

      Fn.add(Fgn.get(0) * -1);
      Fn.add(Fgn.get(1) * -1);
      Fn.add(Fgn.get(2) * -1);

      return Fn;
  }

    private static boolean checkIfMoving(double mewStatic, double Fn, double Fgp)
    {
        Double Ff_max=mewStatic*Fn;
        System.out.println("Ff_max: "+Ff_max);
        if(Ff_max<Fgp)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    private static Vector findFrictionDirection(Vector<Double> Fgp, double fgpLength)
    {
        Vector<Double> vector = new Vector<>();

        double dotProduct = dotProduct(Fgp, Fgp);
        //Multiplying by -1 to make sure it's the opposite direction of the motion
        vector.add(-1*((1/fgpLength)*Fgp.get(0)));
        vector.add(-1*((1/fgpLength)*Fgp.get(1)));
        vector.add(-1*((1/fgpLength)*Fgp.get(2)));
        return vector;

    }

    public static Vector calculateFrictionForce(double fLength, Vector<Double> fDirection)
    {
        Vector<Double> vector = new Vector<>();
        vector.add(fDirection.get(0)*fLength);
        vector.add(fDirection.get(1)*fLength);
        vector.add(fDirection.get(2)*fLength);
        return vector;
    }
    
    public static Vector calculateAcceleration(Object obj, Vector<Double> Fnet)
    {
        Vector<Double> accel=new Vector(3);
        accel.add((1/obj.getMass())*Fnet.get(0));
        accel.add((1/obj.getMass())*Fnet.get(1));
        accel.add((1/obj.getMass())*Fnet.get(2));
        return accel;
    }
    
    public static void eulersMethod(Object obj, Vector<Double> p0, Vector<Double> v0, Vector<Double> fNet, int stepNum, double h, double t0, double tF)
    {
        int step=0;
        
        Vector<Double>a0=calculateAcceleration(obj, fNet);
        System.out.println("          EULERS METHOD");
        System.out.println("*************************************");
        System.out.println("      Initial Conditions");
        System.out.println("P"+step+": "+p0);
        System.out.println("V"+step+": "+v0);
        System.out.println("A"+step+": "+a0);
        
        for(int i=0;i<stepNum;i++)
        {
            step++;
            p0.set(0, (p0.get(0)+(h*v0.get(0))));
            p0.set(1, (p0.get(1)+(h*v0.get(1))));
            p0.set(2, (p0.get(2)+(h*v0.get(2))));
            
            v0.set(0, (v0.get(0)+(h*a0.get(0))));
            v0.set(1, (v0.get(1)+(h*a0.get(1))));
            v0.set(2, (v0.get(2)+(h*a0.get(2))));
            
            calculateAcceleration(obj, fNet);
            
            System.out.println("_______________");
            System.out.println("P"+step+": "+p0);
            System.out.println("V"+step+": "+v0);
            System.out.println("_______________");
        }
        System.out.println("*********************");
    }
}
