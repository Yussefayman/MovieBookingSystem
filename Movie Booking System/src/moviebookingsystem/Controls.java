
package moviebookingsystem;

import static java.awt.SystemColor.window;
import java.io.IOException;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.*;
import javafx.scene.control.Alert;
import javafx.stage.Stage;


public class Controls {
   
    
public  void switchScenes(String dest) throws IOException{
           
          Parent root = FXMLLoader.load(getClass().getResource(dest+".fxml"));
          Scene scene = new Scene(root);
            Main.window.setScene(scene);
            Main.window.setTitle(dest+" Form");
            Main.window.show();
        }
    public static void callAlert(String alertType , String alertTitle , String content){
            Alert alert = null ;
            if(alertType.equals("Error")){
                 alert = new Alert(Alert.AlertType.ERROR);
                alert.setTitle(alertTitle);
                alert.setContentText(content);
                
            }
            else if(alertType.equals("Warning")){
                 alert = new Alert(Alert.AlertType.WARNING);
                alert.setTitle(alertTitle);
                alert.setContentText(content);
                
            }
            else if(alertType.equals("Info")) {
                 alert = new Alert(Alert.AlertType.INFORMATION);
                alert.setTitle(alertTitle);
                alert.setContentText(content);
            }
            alert.showAndWait();
           
        }
}
 

