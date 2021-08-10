
package database;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JOptionPane;


public class mysql {
   
   public static Connection con;
   
   public static void setconnection() throws SQLException{
        
        try {
           Class.forName("com.mysql.jdbc.Driver");
           con= DriverManager.getConnection("jdbc:mysql://localhost:3306/bookingsystem","root","bestprogrammer");
        }
         
         catch (ClassNotFoundException ex) {
            JOptionPane.showMessageDialog(null,ex.getMessage());
        }
    }
    public static boolean runnonquery(String s){
        try{
           mysql.setconnection();
           Statement st=con.createStatement();
            st.execute(s);
            con.close();
            return true;
       }
        catch(SQLException se){
            JOptionPane.showMessageDialog(null,se.getMessage());
            return false;
        }
    }
   public static ResultSet RunOnQuery(String s) throws SQLException{
        try{
           mysql.setconnection();
           Statement stat=con.createStatement();
            ResultSet rs = stat.executeQuery(s);
            
            return rs;
       }
        catch(SQLException se){
           JOptionPane.showMessageDialog(null, se);
           
        }
        con.close();
       return null;
    }
}
