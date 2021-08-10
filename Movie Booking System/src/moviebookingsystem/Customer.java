package moviebookingsystem;

import database.mysql;
import java.sql.ResultSet;
import java.sql.SQLException;


public class Customer {

    private String userName;
    private String password;
    private String email;
    private String id;
    private String phoneNo;
    static int count = 0 ;
  public Customer(){
      
  }
    public Customer(String name, String pass, String email, String phone) {
        this.userName = name;
         this.email = email;
        this.password = pass;
        this.phoneNo = phone;
        count++;
    }
    public void save_data(String name, String email, String password, String phone) throws SQLException{
        
         mysql.runnonquery("insert into customer (  customer_phone, customer_email, customer_username, customer_password) values ('"+phone+"','"+email+"','"+name+"','"+password+"')");
        
    }
    public boolean customerValidation(String Name , String Pass){
        try{   
            ResultSet rs = mysql.RunOnQuery("SELECT * FROM customer");
              while (rs.next()){
                  String username = rs.getString("customer_username");
                  String password = rs.getString("customer_password");
                  if(username.equals(Name) && password.equals(Pass)){
                      //JOptionPane.showMessageDialog(null,"Info INFORMATION : LOGIN LOGGED IN SUCCEFULLY ");
                      return true ;
                  }
              }
           
        }
        catch(Exception e){
            System.out.print(e.getMessage());
           
        }
      return false;
}
    public int customer_id(String name, String password) throws SQLException{
        ResultSet rs = mysql.RunOnQuery("SELECT * FROM customer");
              while (rs.next()){
                  int customer_id= rs.getInt("customer_id");
                 String username = rs.getString("customer_username");
                  String pass = rs.getString("customer_password");
                  
                  if(username.equals(name) && pass.equals(password)){
                    return customer_id;
                  }
    }
              return 0;
    }
}