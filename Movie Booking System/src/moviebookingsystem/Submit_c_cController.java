
package moviebookingsystem;

import java.io.IOException;
import java.net.URL;
import java.util.ResourceBundle;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.FXMLLoader;
import javafx.fxml.Initializable;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.control.ComboBox;
import javafx.stage.Stage;


public class Submit_c_cController implements Initializable {

    @FXML
    private ComboBox<String> city;
    @FXML
    private ComboBox<String> cinema;
    public static String city_;
    public static String cinema_;
    public void submit_c_c(ActionEvent e) throws IOException, Exception {
         city_= city.getValue();
         cinema_= cinema.getValue();
        FXMLLoader loader= new FXMLLoader();
        loader.setLocation(getClass().getResource("films.fxml"));
      
        loader.load();
        FilmsController fc= loader.getController();
        fc.setName(city_, cinema_);
       // fc.pass_c_c(city_, cinema_);
        Parent p= loader.getRoot();
        Stage stage = new Stage ();
        stage.setScene(new Scene(p));
        //Scene scene = new Scene(p);
        stage.showAndWait();
        
        
    }
    
    
    
    @Override
    public void initialize(URL url, ResourceBundle rb) {
       city.getItems().addAll("Cairo","Paris","Alex");
       cinema.getItems().addAll("CarRadio","SuicideIsPainless","AddictWithApen");
    }    
    
}
