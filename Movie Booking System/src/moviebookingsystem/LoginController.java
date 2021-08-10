package moviebookingsystem;
import database.mysql;
import static database.mysql.setconnection;
import java.awt.event.MouseEvent;
import javafx.fxml.FXML;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ResourceBundle;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.*;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.input.MouseButton;
import javax.swing.JOptionPane;

public class LoginController implements Initializable {

    @FXML
    private TextField username;
    @FXML
    private TextField password;
        
    @FXML
     private Button login;
    @FXML
    private Button register;
    public static String username_,pass;
    
    @FXML
    public void actionPerformed(ActionEvent e) throws IOException {
        
         username_ = username.getText();
         pass = password.getText();
        username.setText("");
        password.setText("");
        Customer c = new Customer() ;
        if(checkLogin(username_,pass) && c.customerValidation(username_,pass)){
            Alert nameAlert = new Alert(Alert.AlertType.INFORMATION);
            nameAlert.setTitle("Info: Submission");
                nameAlert.setContentText(" logged in successfully! ");
                nameAlert.showAndWait();
        Controls control = new Controls();
         control.switchScenes("Submit_c_c");
        }
        else if(username_.equals("admin")&&pass.equals("admin")){
             Alert nameAlert = new Alert(Alert.AlertType.INFORMATION);
            nameAlert.setTitle("Info: Submission");
                nameAlert.setContentText(" Logged in Successfully!,Hello Admin ");
                nameAlert.showAndWait();
        Controls control = new Controls();
         control.switchScenes("AdminGui");
        }
    }
    public boolean checkLogin(String name , String pass){
        if (name.length() < 4 || pass.length() < 6) {
            Alert nameAlert = new Alert(Alert.AlertType.WARNING);
                nameAlert.setTitle("WARNING: Adding Failed");
                nameAlert.setContentText(" LOGIN FAILED INVALID PASSWORD ");
                nameAlert.showAndWait();
                return false ;
        }
        return true;
    }
    @FXML
    public void open_signup(ActionEvent e) throws IOException{
         Controls control = new Controls();
         control.switchScenes("RegisterGui");
    }
    @Override
    public void initialize(URL url, ResourceBundle rb) {

    }
  
}