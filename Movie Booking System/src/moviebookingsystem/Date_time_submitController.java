package moviebookingsystem;
import java.io.IOException;
import java.net.URL;
import java.time.LocalDate;
import java.util.ResourceBundle;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.ComboBox;
import javafx.scene.control.DatePicker;



public class Date_time_submitController implements Initializable {
    
     @FXML
     public ComboBox<String> time;
     @FXML
     public ComboBox<String> type;
     @FXML
     public DatePicker date=new DatePicker();
     String s;
    public static LocalDate date_;
    public static String time_,type_;
    @Override
    public void initialize(URL url, ResourceBundle rb) {
       time.getItems().addAll("15:00:00","18:00:00","21:00:00","00:00:00","3:00:00");
       
        

    }  
    @FXML
    public void value_cb(ActionEvent e){
      time_= time.getValue();
      s=  time.getValue();
      if(s.equals("15:00:00")||s.equals("18:00:00")||s.equals("21:00:00"))
       {
           type.getItems().clear();
           type.getItems().addAll("2D");
           date_= date.getValue();
           type_= "2D";
       }
       else if(s.equals("00:00:00")){
           type.getItems().clear();
           type.getItems().addAll("3D");
           date_= date.getValue();
           type_="3D";
       }
       else{
            type.getItems().clear();
           type.getItems().addAll("4D");
           date_= date.getValue();
           type_="4D";
       } 
    }
     @FXML
    public void date_time_submit(ActionEvent e) throws IOException{
         
        Controls control = new Controls();
         control.switchScenes("chair_booking");
         
    }
     
}
