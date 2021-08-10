package moviebookingsystem;

import database.mysql;
import static database.mysql.con;
import java.awt.event.MouseEvent;
import javafx.fxml.FXML;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ResourceBundle;
import java.util.Scanner;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.*;
import javafx.scene.control.Alert.AlertType;
import javax.swing.JOptionPane;

public class RegisterController implements Initializable {
    
 public static String username_,pass;
    Customer c;

    //Function to determine if phoneString Contatin only numbers
    public static boolean isNum(String s) {
        try {
            int num = Integer.parseInt(s);

        } catch (NumberFormatException | NullPointerException npe) {
            return false;
        }
        return true;
    }

    @FXML
    private TextField user;
    @FXML
    private TextField email;
    @FXML
    private TextField password;
    @FXML
    private TextField phone;
   @FXML
     private Button register;

    String username;
    String email_;
    String password_;
    String phone_;
  
  
   ResultSet resultSet=null;
    public boolean valueexist(String user){
        try{
            
          ResultSet rs= mysql.RunOnQuery("select * from customer ");
           while(rs.next()){
               String username= rs.getString("customer_username");
               if(username.equals(user)){
                  
                   return true;
               }
           }
        }
        catch(SQLException se){
          JOptionPane.showMessageDialog(null, se);
        }
        return false;
    }
    
    
    @FXML
    public boolean check(){
        try{
        boolean u_nconfirmed = true;
        boolean emailconfirmed = true;
        boolean passconfirmed = true;
        boolean phoneconfirmed = true;
        boolean iDconfirmed = true;
        boolean unAlreadyExist = false;
        if (username.length() < 4) {
            u_nconfirmed = false;
            Alert nameAlert = new Alert(AlertType.WARNING);
            nameAlert.setTitle("WARNING: Registration Failed");
            nameAlert.setContentText("Username is too short ");
            nameAlert.showAndWait();
        }
        if(valueexist(username)){
            unAlreadyExist = false;
             Alert nameAlert = new Alert(AlertType.WARNING);
                nameAlert.setTitle("WARNING: Registration Failed");
                nameAlert.setContentText("Username is already exist :( ");
                nameAlert.showAndWait();
        }
         if (!email.getText().contains("@gmail.com") && !email.getText().contains("@hotmail.com") && !email.getText().contains("@yahoo.com")) {
            emailconfirmed = false;
            Alert emailAlert = new Alert(AlertType.WARNING);
            emailAlert.setTitle("WARNING: Registration Failed");
            emailAlert.setContentText("Email is invalid ");
            emailAlert.showAndWait();
        }
        if (password.getText().length() < 8) {
            passconfirmed = false;
            Alert passAlert = new Alert(AlertType.WARNING);
            passAlert.setTitle("WARNING: Registration Failed");
            passAlert.setContentText("Password is too short ");
            passAlert.showAndWait();
        }
         if (!isNum(phone.getText()) || phone.getText().length() < 10) {
            phoneconfirmed = false;
            Alert pAlert = new Alert(AlertType.WARNING);
            pAlert.setTitle("WARNING: Registration Failed");
            pAlert.setContentText("Please enter correct phone number !");
            pAlert.showAndWait();
        }
        
         if (emailconfirmed && passconfirmed && u_nconfirmed && phoneconfirmed && iDconfirmed && !unAlreadyExist) {
            Alert rAlert = new Alert(AlertType.INFORMATION);
            rAlert.setTitle("Registerition Info ");
            rAlert.setContentText("You've Registered Succefully !");
            rAlert.showAndWait();
            return true;
         }
        }
        catch(Exception se){
            JOptionPane.showMessageDialog(null, se);
        }
        return false;
    }
    
    
    @FXML
     public void handleButtonAction(ActionEvent me) throws SQLException, IOException{
         
         if( Registration())
             
         { Controls control = new Controls();
           control.switchScenes("submit_c_c");
         }
     }
    
    public boolean Registration() throws SQLException{
      
          username=user.getText();
          email_=email.getText();
          password_=password.getText();
          phone_=phone.getText();   
         if (check()){
             username_=user.getText();
             pass= password.getText();
         Customer c= new Customer(username, email_, password_, phone_);
         c.save_data(username, email_, password_, phone_);
          return true;
         }
         return false;
      }
     
    @Override
    public void initialize(URL url, ResourceBundle rb) {

    }

}
