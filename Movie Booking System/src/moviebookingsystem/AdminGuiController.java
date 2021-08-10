
package moviebookingsystem;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.sql.SQLException;
import java.util.ResourceBundle;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.Node;
import javafx.scene.control.Alert;
import javafx.scene.control.ComboBox;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javax.imageio.ImageIO;


public class AdminGuiController implements Initializable {

    @FXML
    private TextField name;
    @FXML
    private TextField genre;
    @FXML
    private TextArea description;
    @FXML
    private TextField release_date;
     @FXML
    private ComboBox<String> city;
    @FXML
    private ComboBox<String> cinema;
    @FXML
    private ImageView imageView;
    @FXML
    private TextField hall_id;
     
    private String name_,releasedate_,genre_,description_,poster_,city_,cinema_,hallid;
    private FileChooser fc;
    private File filepath;
    private char [] imageName;
    String path;
    
    @FXML
    public void chooseimage(ActionEvent e) throws IOException{
      Stage stage=(Stage)((Node)e.getSource()).getScene().getWindow();
      fc= new FileChooser();
      fc.setTitle("open image");
      this.filepath=fc.showOpenDialog(stage);
      
      String userdirect= System.getProperty("user.home");
      File userDirectory= new File(userdirect);
      if(!userDirectory.canRead())
          userDirectory=new File("D:/Pictures");
          fc.setInitialDirectory(userDirectory);
         
              BufferedImage img = ImageIO.read(filepath);
               path= filepath.getAbsolutePath();
              Image image= SwingFXUtils.toFXImage(img, null);
              imageView.setImage(image);
              
     
    }
    
    
     @FXML
     public void savedata(ActionEvent e) throws SQLException{
         imageName = new char [10];
         
         name_=name.getText();
         releasedate_=release_date.getText();
         genre_=genre.getText();
         description_=description.getText();
         //poster_=s;
          String city_= city.getValue();
           String cinema_= cinema.getValue();
         hallid=hall_id.getText();
         int counter = 0 ;
         int pathlength = path.length();
         for (int i = pathlength -1  ; i > 0 ; i--){
             if(path.charAt(i)=='\\') break;
             
             else {
                  imageName[counter]=path.charAt(i);
                  counter++;
             }
           
         }
         int count = 0 ;
         
         char [] x = new char [10];
         for(int i =imageName.length -1 ; i >=0 ; i--){
             x[count]=imageName[i];
             count++;
         }
         path = new String(x);
         System.out.print(path);
         
         movie mo= new movie(name_,releasedate_,genre_,description_,path);
         cinema c= new cinema(cinema_, city_);
         hall h= new hall(Integer.valueOf(hallid));
         c.save_cinema(cinema_, city_);
         if(c.save_data(cinema_, city_, name_, releasedate_ ,genre_ , description_, path, Integer.valueOf(hallid))){
            
         Controls.callAlert("Info", "CONFIRM INFORMATION", "The movie has been submitted successfully");
         }
         else{
            Alert nameAlert = new Alert(Alert.AlertType.WARNING);
                nameAlert.setTitle("WARNING: Movie Saving Failed");
                nameAlert.setContentText("Movie isn't saved :( ");
                nameAlert.showAndWait();
         }
     }
    
    @Override
    public void initialize(URL url, ResourceBundle rb) {
        city.getItems().addAll("Cairo","Paris","Alex");
       cinema.getItems().addAll("CarRadio","SuicideIsPainless","AddictWithApen");
    }    
    
}
